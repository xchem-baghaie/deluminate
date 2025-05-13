from typing import List, Tuple
from itertools import chain

from rdkit import Chem

from tools.parallel_map import parallel_map_generic
from transforms.rdkit_mols import compute_smiles_from_mol

from configs.stereochemistry import StereoInfo, BRIDGED_SMARTS_STRING

bridged_smarts = Chem.MolFromSmarts(BRIDGED_SMARTS_STRING)


def calculate_num_stereoisomers(mol_stereoname: Tuple) -> int:
    """
    Given a tuple: (mol (type Chem.rdchem.Mol), stereochemistryname (type string)), calculates the number of possible stereoisomers which can physically exist.
    If the stereochemistryname indicates that the number of isomers is double what can naturally be captured by the SMILES string, it returns the more appropriate (doubled) value.
    """
    mol = mol_stereoname[0]
    stereoname = mol_stereoname[1]
    if compute_smiles_from_mol(mol) == "":
        return 0
    else:
        si = StereoInfo(mol)
        stereoisomer_count = (
            si.get_stereoisomer_count_only_physically_plausible_varying_only_undefined_stereocenters()
        )
        if stereoname in [
            "Relative Cis",
            "Relative Trans",
            "Syn",
            "Anti",
            "Regioisomeric Mixture",
        ]:
            if not si.is_meso():
                return (
                    stereoisomer_count * 2
                )  # There are double the amount of stereoisomers present, which cannot be captured by the SMILES string
            return stereoisomer_count  # net result of doubling amount of stereoisomers present, then halving because it is meso
        elif stereoname in ["Endo", "Exo", "Endo/Exo"]:
            # Determine whether stereochemistry was specified at the bridgehead carbons...if it was, then we need to double the return value. If it was not,
            # we do not have to do this, because RDKit should have enumerated both conformations anyways.
            if si.has_undefined_bridgehead(bridged_smarts):
                return stereoisomer_count
            return stereoisomer_count * 2
        else:
            if not si.is_meso():
                return stereoisomer_count  # BB does not fall into any special case, based on its stereoname
            return int(
                stereoisomer_count / 2
            )  # BB has half the expected number of stereoisomers because it is meso


def parallel_map_calculate_num_stereoisomers(
    mol_array: List[Chem.rdchem.Mol], stereoname_array: List[str]
) -> List[int]:
    """
    Parallel map implementation of calculate_num_stereoisomers
    """

    if len(mol_array) != len(stereoname_array):
        raise Exception(
            f"Arrays of differing lengths ({len(mol_array)} and {len(stereoname_array)}) were passed into parallel_map_calculate_num_stereoisomers"
        )

    mol_stereoname_array = [
        (mol_array[i], stereoname_array[i]) for i in range(len(mol_array))
    ]

    return list(
        chain.from_iterable(
            parallel_map_generic(mol_stereoname_array, calculate_num_stereoisomers)
        )
    )
