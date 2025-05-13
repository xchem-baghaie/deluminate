import itertools
from copy import copy

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)
from rdkit.Chem.MolKey.InchiInfo import InchiInfo
from rdkit.Chem.rdmolops import FindPotentialStereo
from rdkit.Chem.rdchem import StereoType, ChiralType

from configs import Config

MAX_N_STEREOISOMERS_ALLOWED = 4
ALLENELIKE_SMARTS = "[AD2](=A)=A"
BRIDGED_SMARTS_STRING = "[$(A~1~2~*~A(~*~1)~*~2),$(A~1~2~*~A(~*~1)~*~*~2),$(A~1~2~*~A(~*~1)~*~*~*~2),$(A~1~2~*~A(~*~1)~*~*~*~*~2),$(A~1~2~*~A(~*~*~1)~*~*~2),$(A~1~2~*~A(~*~*~1)~*~*~*~2),$(A~1~2~*~A(~*~*~1)~*~*~*~*~2),$(A~1~2~*~A(~*~*~*~1)~*~*~*~2),$(A~1~2~*~A(~*~*~*~1)~*~*~*~*~2),$(A~1~2~*~A(~*~*~*~*~1)~*~*~*~*~2),$(A~1~2~*~*~A(~*~*~1)~*~*~2),$(A~1~2~*~*~A(~*~*~1)~*~*~*~2),$(A~1~2~*~*~A(~*~*~1)~*~*~*~*~2),$(A~1~2~*~*~A(~*~*~*~1)~*~*~*~2),$(A~1~2~*~*~A(~*~*~*~1)~*~*~*~*~2),$(A~1~2~*~*~A(~*~*~*~*~1)~*~*~*~*~2),$(A~1~2~*~*~*~A(~*~*~*~1)~*~*~*~2),$(A~1~2~*~*~*~A(~*~*~*~1)~*~*~*~*~2),$(A~1~2~*~*~*~A(~*~*~*~*~1)~*~*~*~*~2),$(A~1~2~*~*~*~*~A(~*~*~*~*~1)~*~*~*~*~2)]"
RACEMIC_AXIAL_SMARTS_QUERIES = [
    "[CD2](=[$([C!H2]),$([N!H1])])=[$([C!H2]),$([N!H1])]",  # Allene, Carbodiimide
    "[C!H2]1[CH2][CH0]2([CH2][C!H2][CH2]2)[CH2]1",  # spiro[3.3]heptane
    "[C!H2]1[CH2][CH0]2([CH2][CH2][C!H2][CH2][CH2]2)[CH2]1",  # spiro[3.4]nonane
    "[C!H2]1[CH2][CH2][CH0]2([CH2][CH2][C!H2][CH2][CH2]2)[CH2][CH2]1",  # spiro[4.4]undecane
    "[C!H2]1[CH2][CH2][CH0]2([CH2][CH2][CH2][C!H2][CH2][CH2][CH2]2)[CH2][CH2]1",  # spiro[4.5]triundecane
    "[C!H2]1[CH2][CH2][CH2][CH0]2([CH2][CH2][CH2][C!H2][CH2][CH2][CH2]2)[CH2][CH2][CH2]1",  # spiro[5.5]pentadecane
    "[C!H2]1[CH2][CH0]2([CH2][CH2][CH2][C!H2][CH2][CH2][CH2]2)[CH2]1",  # spiro[3.5]undecane
]

STEREONAMES = [
    "Absolute",
    "Achiral",
    "Syn",
    "Anti",
    "Diastereomeric Mixture",
    "E",
    "E/Z",
    "Endo",
    "Endo/Exo",
    "Exo",
    "Multiple Stereochemical Relationships",
    "Racemic",
    "Racemic (Axial)",
    "Relative Cis",
    "Relative Trans",
    "Z",
]

bridged_smarts = Chem.MolFromSmarts(BRIDGED_SMARTS_STRING)


class StereoInfo(Config):
    """
    A class that describes inherent stereochemical properties of molecules.
    All class properties are lists of 2-tuples, where each element in the list corresponds
    to a single stereogenic atom or bond. The first element in each tuple is the id of
    the stereogenic atom or bond, and the second element in the tuple is the stereochemical
    descriptor of that atom or bond.

    Properties:
    - AtomStereocenters_All: List[Tuple]
        Any atoms that have tetrahedral chirality
    - AtomStereocenters_Para: List[Tuple]
        Any atoms that have tetrahedral, para stereochemistry (e.g. 1,4-dimethyl cyclohexane)
    - AtomStereocenters_Not_Para: List[Tuple]
        Any atoms that have tetrahedral chirality that are not part of a para group
    - DefinedAtomStereocenters_All: List[Tuple]
        Any atoms that have tetrahedral chirality that are defined
    - DefinedAtomStereocenters_Para: List[Tuple]
        Any atoms that have tetrahedral, para stereochemistry that are defined
    - DefinedAtomStereocenters_Not_Para: List[Tuple]
        Any atoms that have tetrahedral chirality that are not part of a para group and are defined
    - UndefinedAtomStereocenters_All: List[Tuple]
        Any atoms that have tetrahedral chirality that are undefined
    - UndefinedAtomStereocenters_Para: List[Tuple]
        Any atoms that have tetrahedral, para stereochemistry that are undefined
    - UndefinedAtomStereocenters_Not_Para: List[Tuple]
        Any atoms that have tetrahedral chirality that are not part of a para group and are undefined
    - BondStereocenters: List[Tuple]
        Any bonds that are stereogenic
    - DefinedBondStereocenters: List[Tuple]
        Any bonds that are stereogenic and defined
    - UndefinedBondStereocenters: List[Tuple]
        Any bonds that are stereogenic and undefined

    Methods:
    - is_meso
        Args:
        Returns: Bool (whether or not the molecule is meso)
        Raises:
    - contains_racemic_axial_substructure
        Args: spiro_3_3_heptane_smarts: str, spiro_4_4_undecane_smarts: str, disubstituted_allene_smarts: str
        Returns Bool (whether or not the molecule matches either of the SMARTS string queries)
        Raises:
    - has_bridge
        Args: bridged_smarts: str
        Returns: Bool (whether or not the molecule has a bridged ring system)
        Raises:
    - has_undefined_bridgehead
        Args: bridged_smarts: str
        Returns: Bool (whether or not the molecule has undefined stereochemistry at a bridgehead atom)
        Raises:
    - has_bridge_but_is_neither_endo_nor_exo
        Args: bridged_smarts: str
        Returns: Bool (whether or not the molecule has a bridged ring system, but is neither endo nor exo)
        Raises:
    - has_para_stereo
        Args:
        Returns: Bool (whether or not the molecule has para stereochemistry)
        Raises:
    - get_stereoisomer_count_only_physically_plausible_varying_only_undefined_stereocenters
        Args:
        Returns: Int (name describes return value)
        Raises:
    - get_stereoisomer_count_only_physically_plausible_varying_all_stereocenters
        Args:
        Returns: Int (name describes return value)
        Raises:
    - get_stereoisomer_count_no_physical_plausibility_check_varying_only_undefined_stereocenters
        Args:
        Returns: Int (name describes return value)
        Raises:
    - get_stereoisomer_count_no_physical_plausibility_check_varying_all_stereocenters
        Args:
        Returns: Int (name describes return value)
        Raises:
    - contains_physically_restrictive_group
        Args:
        Returns: Bool (whether the molecule has a group that restricts the orientation of some stereogenic atoms)
        Raises:
    - db_directionality
        Args:
        Returns: Str (If the molecule has exactly 1 isomerizable double bond, the description of the double bond directionality. Otherwise, it returns 'Multiple' or None)
        Raises:
    - all_stereocenters_in_same_ring
        Args:
        Returns: Bool: (Whether or not all stereocenters are in the same ring)
        Raises:
    - all_stereocenters_acyclic
        Args:
        Returns: Bool: (Whether or not all stereocenters are on atoms which are acyclic)
        Raises:
    - all_stereocenters_cyclic
        Args:
        Returns: Bool: (Whether or not all stereocenters are on atoms which are cyclic)
        Raises:
    - flag_stereoname
        Args: stereoname: str
        Returns: Bool (Whether or not the stereochemistry name provided is consistent with the structure as drawn)
        Raises:
    """

    def __init__(self, mol):
        self.mol = mol

        """ Atom Stereocenters (tetrahedral chirality) """
        self.AtomStereocenters_All = Chem.FindMolChiralCenters(
            self.mol, force=True, includeUnassigned=True, useLegacyImplementation=False
        )
        # Define all undefined tetrahedral stereocenters to id which ones are para stereocenters
        tmp_mol = copy(self.mol)
        for x in self.AtomStereocenters_All:
            if x[1] == "?":
                tmp_mol.GetAtomWithIdx(x[0]).SetChiralTag(
                    Chem.ChiralType.CHI_TETRAHEDRAL_CW
                )
        tmp_stereo = Chem.FindMolChiralCenters(
            tmp_mol, force=True, includeUnassigned=True, useLegacyImplementation=False
        )
        tmp_stereo_dct = {
            self.AtomStereocenters_All[i]: tmp_stereo[i]
            for i in range(len(self.AtomStereocenters_All))
        }
        self.AtomStereocenters_Para = [
            x
            for x in self.AtomStereocenters_All
            if (tmp_stereo_dct[x][1] in ["r", "s"] or x[1] in ["r", "s"])
        ]
        self.AtomStereocenters_Not_Para = [
            x
            for x in self.AtomStereocenters_All
            if (tmp_stereo_dct[x][1] not in ["r", "s"] and x[1] not in ["r", "s"])
        ]
        assert len(self.AtomStereocenters_All) == len(
            self.AtomStereocenters_Para
        ) + len(self.AtomStereocenters_Not_Para)

        self.DefinedAtomStereocenters_All = [
            x
            for x in self.AtomStereocenters_All
            if x[1] in ["R", "S", "r", "s", "Tet_CW", "Tet_CCW"]
        ]
        self.DefinedAtomStereocenters_Para = [
            x for x in self.AtomStereocenters_Para if x[1] != "?"
        ]
        self.DefinedAtomStereocenters_Not_Para = [
            x for x in self.AtomStereocenters_Not_Para if x[1] != "?"
        ]
        assert len(self.DefinedAtomStereocenters_All) == len(
            self.DefinedAtomStereocenters_Para
        ) + len(self.DefinedAtomStereocenters_Not_Para)

        self.UndefinedAtomStereocenters_All = [
            x for x in self.AtomStereocenters_All if x[1] == "?"
        ]
        self.UndefinedAtomStereocenters_Para = [
            x for x in self.AtomStereocenters_Para if x[1] == "?"
        ]
        self.UndefinedAtomStereocenters_Not_Para = [
            x for x in self.AtomStereocenters_Not_Para if x[1] == "?"
        ]
        assert len(self.UndefinedAtomStereocenters_All) == len(
            self.UndefinedAtomStereocenters_Para
        ) + len(self.UndefinedAtomStereocenters_Not_Para)

        assert len(self.AtomStereocenters_All) == len(
            self.DefinedAtomStereocenters_All
        ) + len(self.UndefinedAtomStereocenters_All)
        assert len(self.AtomStereocenters_Para) == len(
            self.DefinedAtomStereocenters_Para
        ) + len(self.UndefinedAtomStereocenters_Para)
        assert len(self.AtomStereocenters_Not_Para) == len(
            self.DefinedAtomStereocenters_Not_Para
        ) + len(self.UndefinedAtomStereocenters_Not_Para)

        """ Bond Stereocenters """
        db_dct = {
            Chem.rdchem.BondStereo.STEREOCIS: "Z",
            Chem.BondStereo.STEREOZ: "Z",
            Chem.rdchem.BondStereo.STEREOTRANS: "E",
            Chem.BondStereo.STEREOE: "E",
            Chem.BondStereo.STEREOANY: "Undefined",
            Chem.BondStereo.STEREONONE: "Undefined",
        }
        self.BondStereocenters = [
            (si.centeredOn, db_dct[self.mol.GetBondWithIdx(si.centeredOn).GetStereo()])
            for si in FindPotentialStereo(self.mol)
            if si.type == StereoType.Bond_Double
        ]
        self.DefinedBondStereocenters = [
            x for x in self.BondStereocenters if x[1] in ["E", "Z"]
        ]
        self.UndefinedBondStereocenters = [
            x for x in self.BondStereocenters if x[1] == "Undefined"
        ]
        assert len(self.BondStereocenters) == len(self.DefinedBondStereocenters) + len(
            self.UndefinedBondStereocenters
        )

    def is_meso(self):
        return InchiInfo(Chem.MolToInchi(self.mol)).get_sp3_stereo()["main"][
            "non-isotopic"
        ][2]

    def contains_racemic_axial_substructure(
        self, smarts_query_list=RACEMIC_AXIAL_SMARTS_QUERIES
    ):
        flag = False  # init
        for query in smarts_query_list:
            smarts = Chem.MolFromSmarts(query)
            if self.mol.HasSubstructMatch(smarts):
                flag = True
                break
        return flag

    def has_bridge(self, bridged_smarts: str):
        return self.mol.HasSubstructMatch(bridged_smarts)

    def has_undefined_bridgehead(self, bridged_smarts: str):
        hasundefinedbridgehead = False
        atoms = self.mol.GetSubstructMatches(bridged_smarts)
        atomidxs = []
        for atomno in range(len(atoms)):
            atomidxs.append(atoms[atomno][0])
        for atomidx in atomidxs:
            if (
                self.mol.GetAtomWithIdx(atomidx).GetChiralTag()
                == ChiralType.CHI_UNSPECIFIED
            ):
                hasundefinedbridgehead = True
        return hasundefinedbridgehead

    def has_bridge_but_is_neither_endo_nor_exo(self, bridged_smarts: str):
        if not self.has_bridge(bridged_smarts):
            return False
        else:
            stereocentersnotatbridges = [
                x
                for x in self.AtomStereocenters_All
                if x not in self.mol.GetSubstructMatches(bridged_smarts)
            ]
            atomidxs = [x[0] for x in stereocentersnotatbridges]
            parities = [x[1] for x in stereocentersnotatbridges]
            inlessthan3rings_ls = []
            ri = self.mol.GetRingInfo()
            for idx in atomidxs:
                inlessthan3rings_ls.append(ri.NumAtomRings(idx) < 3)
            if parities.count("?") > 1 and all(inlessthan3rings_ls):
                return True
            else:
                return False

    def has_para_stereo(self):
        return len(self.AtomStereocenters_Para) > 0

    def get_stereoisomer_count_only_physically_plausible_varying_only_undefined_stereocenters(
        self,
    ):
        return len(
            tuple(
                EnumerateStereoisomers(
                    self.mol,
                    options=StereoEnumerationOptions(
                        tryEmbedding=True, unique=True, onlyUnassigned=True
                    ),
                )
            )
        )

    def get_stereoisomer_count_only_physically_plausible_varying_all_stereocenters(
        self,
    ):
        return len(
            tuple(
                EnumerateStereoisomers(
                    self.mol,
                    options=StereoEnumerationOptions(
                        tryEmbedding=True, unique=True, onlyUnassigned=False
                    ),
                )
            )
        )

    def get_stereoisomer_count_no_physical_plausibility_check_varying_only_undefined_stereocenters(
        self,
    ):
        return len(
            tuple(
                EnumerateStereoisomers(
                    self.mol,
                    options=StereoEnumerationOptions(
                        tryEmbedding=False, unique=True, onlyUnassigned=True
                    ),
                )
            )
        )

    def get_stereoisomer_count_no_physical_plausibility_check_varying_all_stereocenters(
        self,
    ):
        return len(
            tuple(
                EnumerateStereoisomers(
                    self.mol,
                    options=StereoEnumerationOptions(
                        tryEmbedding=False, unique=True, onlyUnassigned=False
                    ),
                )
            )
        )

    def contains_physically_restrictive_group(self):
        return (
            self.get_stereoisomer_count_no_physical_plausibility_check_varying_all_stereocenters()
            != self.get_stereoisomer_count_only_physically_plausible_varying_all_stereocenters()
        )

    def db_directionality(self):
        if len(self.BondStereocenters) == 0:
            return None
        elif len(self.BondStereocenters) >= 2:
            return "Multiple"
        elif self.mol.HasSubstructMatch(Chem.MolFromSmarts("[AD2](=A)=A")):
            return None
        else:
            return self.BondStereocenters[0][1]

    def all_stereocenters_in_same_ring(self):
        tet_stereo_para_and_not_para_atom_idxs = [
            x[0] for x in self.AtomStereocenters_All
        ]
        ri = self.mol.GetRingInfo()
        bool_list = []
        for combo in itertools.combinations(tet_stereo_para_and_not_para_atom_idxs, 2):
            bool_list.append(ri.AreAtomsInSameRing(combo[0], combo[1]))
        return all(bool_list)

    def all_stereocenters_acyclic(self):
        tet_stereo_para_and_not_para_atom_idxs = [
            x[0] for x in self.AtomStereocenters_All
        ]
        ri = self.mol.GetRingInfo()
        bool_list = []
        for atomidx in tet_stereo_para_and_not_para_atom_idxs:
            bool_list.append(ri.NumAtomRings(atomidx) == 0)
        return all(bool_list)

    def all_stereocenters_cyclic(self):
        tet_stereo_para_and_not_para_atom_idxs = [
            x[0] for x in self.AtomStereocenters_All
        ]
        ri = self.mol.GetRingInfo()
        bool_list = []
        for atomidx in tet_stereo_para_and_not_para_atom_idxs:
            bool_list.append(ri.NumAtomRings(atomidx) > 0)
        return all(bool_list)

    def flag_stereoname(self, stereoname: str):
        n_allenelike_centers = len(
            self.mol.GetSubstructMatches(Chem.MolFromSmarts(ALLENELIKE_SMARTS))
        )
        total_stereo = (
            len(self.AtomStereocenters_All)
            + len(self.BondStereocenters)
            - n_allenelike_centers
        )
        has_tet_stereo = len(self.AtomStereocenters_All) > 0
        has_db_stereo = len(self.BondStereocenters) > 0
        mixed_tet_and_db_stereo = (
            len(self.BondStereocenters) > 0 and len(self.AtomStereocenters_All) > 0
        )
        num_tet_stereo_para_and_not_para = len(self.AtomStereocenters_All)
        num_def_tet_stereo_para_and_not_para = len(self.DefinedAtomStereocenters_All)
        num_undef_tet_stereo_para_and_not_para = len(
            self.UndefinedAtomStereocenters_All
        )
        num_db_stereo = len(self.BondStereocenters) - n_allenelike_centers
        num_def_db_stereo = len(self.DefinedBondStereocenters)
        num_undef_db_stereo = (
            len(self.UndefinedBondStereocenters) - n_allenelike_centers
        )
        if num_db_stereo == 0:  # Reset in case allenelike substructures were present
            has_db_stereo = False
            mixed_tet_and_db_stereo = False
        db_directionality = self.db_directionality()
        has_para = len(self.AtomStereocenters_Para) > 0
        contains_bridge = self.has_bridge(bridged_smarts)
        all_stereocenters_in_same_ring = self.all_stereocenters_in_same_ring
        has_bridge_but_is_neither_endo_nor_exo = (
            self.has_bridge_but_is_neither_endo_nor_exo(bridged_smarts)
        )
        all_stereocenters_acyclic = self.all_stereocenters_acyclic()
        all_stereocenters_cyclic = self.all_stereocenters_cyclic()
        has_racemic_axial_substructure = self.contains_racemic_axial_substructure()
        if stereoname == "":
            return True
        elif stereoname not in STEREONAMES:
            return True
        elif stereoname == "Absolute":
            return not (
                total_stereo >= 1
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and not has_para
            )
        elif stereoname == "Achiral":
            return not (total_stereo == 0)
        elif stereoname == "Syn":
            return not (
                total_stereo == 2
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and not has_para
                and all_stereocenters_acyclic
            )
        elif stereoname == "Anti":
            return not (
                total_stereo == 2
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and not has_para
                and all_stereocenters_acyclic
            )
        elif stereoname == "Diastereomeric Mixture":
            return not (
                total_stereo >= 2
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para < total_stereo
                and num_undef_tet_stereo_para_and_not_para >= 1
                and (not contains_bridge or has_bridge_but_is_neither_endo_nor_exo)
                and not has_racemic_axial_substructure
            )
        elif stereoname == "E":
            return not (
                total_stereo == 1
                and not has_tet_stereo
                and has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_db_stereo == 1
                and num_def_db_stereo == 1
                and num_undef_db_stereo == 0
                and db_directionality == "E"
                and not has_para
            )
        elif stereoname == "E/Z":
            return not (
                total_stereo == 1
                and not has_tet_stereo
                and has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_db_stereo == 1
                and num_def_db_stereo == 0
                and num_undef_db_stereo == 1
                and db_directionality == "Undefined"
                and not has_para
            )
        elif stereoname == "Endo":
            return not (
                total_stereo >= 3
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and not has_para
                and contains_bridge
                and not has_bridge_but_is_neither_endo_nor_exo
            )
        elif stereoname == "Endo/Exo":
            return not (
                total_stereo >= 3
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para < total_stereo
                and num_undef_tet_stereo_para_and_not_para >= 1
                and not has_para
                and contains_bridge
                and not has_bridge_but_is_neither_endo_nor_exo
            )
        elif stereoname == "Exo":
            return not (
                total_stereo >= 3
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and not has_para
                and contains_bridge
                and not has_bridge_but_is_neither_endo_nor_exo
            )
        elif stereoname == "Multiple Stereochemical Relationships":
            return not (
                (
                    # All tetrahedral case
                    total_stereo >= 3
                    and has_tet_stereo
                    and not has_db_stereo
                    and not mixed_tet_and_db_stereo
                    and num_def_tet_stereo_para_and_not_para == total_stereo
                    and not contains_bridge
                )
                or (
                    # All double bond case
                    not has_tet_stereo
                    and has_db_stereo
                    and not mixed_tet_and_db_stereo
                    and num_db_stereo >= 2
                )
                or (
                    # Mixed tetrahedral and double bond case
                    has_tet_stereo
                    and has_db_stereo
                    and mixed_tet_and_db_stereo
                )
            )
        elif stereoname == "Racemic":
            return not (
                total_stereo == 1
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == 1
                and num_def_tet_stereo_para_and_not_para == 0
                and num_undef_tet_stereo_para_and_not_para == 1
                and not has_para
            )
        elif stereoname == "Racemic (Axial)":
            return not (
                (
                    (
                        total_stereo == 3
                        and has_tet_stereo
                        and not has_db_stereo
                        and not mixed_tet_and_db_stereo
                        and num_tet_stereo_para_and_not_para == 3
                        and num_def_tet_stereo_para_and_not_para == 0
                        and num_undef_tet_stereo_para_and_not_para == 3
                        and not has_para
                    )
                    or (
                        total_stereo == 1
                        and not has_tet_stereo
                        and has_db_stereo
                        and num_def_db_stereo == 0
                        and num_undef_db_stereo == 1
                    )
                )
                and has_racemic_axial_substructure
            )
        elif stereoname == "Regioisomeric Mixture":
            return False
        elif stereoname == "Relative Cis":
            return not (
                total_stereo == 2
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and all_stereocenters_in_same_ring
                and all_stereocenters_cyclic
            )
        elif stereoname == "Relative Trans":
            return not (
                total_stereo == 2
                and has_tet_stereo
                and not has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_tet_stereo_para_and_not_para == total_stereo
                and num_def_tet_stereo_para_and_not_para == total_stereo
                and num_undef_tet_stereo_para_and_not_para == 0
                and all_stereocenters_in_same_ring
                and all_stereocenters_cyclic
            )
        elif stereoname == "Z":
            return not (
                total_stereo == 1
                and not has_tet_stereo
                and has_db_stereo
                and not mixed_tet_and_db_stereo
                and num_db_stereo == 1
                and num_def_db_stereo == 1
                and num_undef_db_stereo == 0
                and db_directionality == "Z"
                and not has_para
            )
        else:
            return True
