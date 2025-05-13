from copy import deepcopy
from dataclasses import dataclass
from typing import Tuple, List, Dict, Optional, Union

import math
import numpy as np
import pandas as pd
import pandas.io.formats.excel
import pyarrow as pa
from pyarrow import dataset as pda
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

from configs import Config
from configs.xviewer import XViewer_Connection

@dataclass
class FunctionalGroupConfig(Config):
    """
    Attributes:
    ___________

    fg_name:
        -str, the name of the functional group
    smarts:
        -str, the SMARTS string used to query for the functional group
    example_smiles:
        -str, the SMILES string of an example structure which contains this functional group
    compound_proposal_rule:
        -str, one of "Ignore", "Warn", or "Alert". Can be used to filter out undesirable functional groups (e.g. "Alert"), or flag them (e.g. "Warn")
    new_bb_purchase_rule:
        -str, one of "Ignore" or "Alert". Any "Alert" functional groups would not be considered for purchase.
    xviewer_connection:
        -XViewer_Connection, the connection to the XViewer database
            """

    def __init__(
        self,
        fg_name: str,
        smarts: str,
        example_smiles: str,
        compound_proposal_rule: str,
        new_bb_purchase_rule: str,
        xviewer_connection: XViewer_Connection = None,
    ):
        self.fg_name = fg_name
        self.smarts = smarts
        self.example_smiles = example_smiles
        self.compound_proposal_rule = compound_proposal_rule
        self.new_bb_purchase_rule = new_bb_purchase_rule
        self.xviewer_connection = xviewer_connection


@dataclass
class ReactionConfig(Config):
    """
    Attributes:
    ___________

    rxn_name:
        -str, the name of the reaction (e.g. "Reductive Amination")
    smarts:
        -str, the SMARTS representation of the reaction
    allowed_fgs_reactant1:
        -List[FunctionalGroupConfig], the list of functional groups in the first reactant which can participate in the reaction
    forbidden_fgs_reactant1:
        -List[FunctionalGroupConfig], the list of functional groups in the first reactant which are incompatible with the reaction
    example_reactant1:
        -str, the SMILES string of an example molecule which may be used as the first reactant in this reaction
    intermediate_isolated:
        -bool, whether any intermediate is isolated over the course of the overall transformation represented by the SMARTS string
    associated_validation_assays_reactant1 (optional):
        -str, the name of the primary assay which is used to validate the first reactant in this reaction
    associated_validation_assays_reactant2 (optional):
        -str, the name of the primary assay which is used to validate the second reactant in this reaction
    allowed_fgs_reactant2 (optional):
        -List[FunctionalGroupConfig], the list of functional groups in the second reactant in this reaction
    forbidden_fgs_reactant2 (optional):
        -List[FunctionalGroupConfig], the list of functional groups in the second reactant which are incompatible with the reaction
    product_fgs:
        -List[FunctionalGroupConfig], the list of functional groups that are present in the product of the reaction (not all will necessarily be present, but at least one will)
    example_reactant2 (optional):
        -str, the SMILES string of an example molecule which may be used as the second reactant in this reaction
    example_product (optional):
        -str, the SMILES string of an example molecule which would be the product of the reaction
    manually_added_fgs_reactant2 (optional):
        -List[FunctionalGroupConfig], any functional groups the user would like to include in the second reactant
    manually_excluded_fgs_reactant2 (optional):
        -List[FunctionalGroupConfig], any functional groups the user would like to exclude from the second reactant
    branch_step (optional):
        -Tuple, of the form (reaction branch idx, step) where both elements of the tuple are integers.
    num_reactants:
        -int, the number of reactants that participate in the reaction. This is calculated from the SMARTS string, as the number of period characters + 1. This may be overwritten for special cases.
    xviewer_connection:
        -XViewer_Connection, the connection to the XViewer database
    """

    def __init__(
        self,
        rxn_name: str,
        smarts: str,
        allowed_fgs_reactant1: List[FunctionalGroupConfig],
        forbidden_fgs_reactant1: List[FunctionalGroupConfig],
        example_reactant1: str,
        intermediate_isolated: bool,
        associated_validation_assays_reactant_1: Optional[List[str]] = None,
        associated_validation_assays_reactant_2: Optional[List[str]] = None,
        allowed_fgs_reactant2: Optional[List[FunctionalGroupConfig]] = [],
        forbidden_fgs_reactant2: Optional[List[FunctionalGroupConfig]] = [],
        product_fgs: Optional[List[FunctionalGroupConfig]] = [],
        example_reactant2: Optional[str] = None,
        example_product: Optional[str] = None,
        manually_added_fgs_reactant2: Optional[List[FunctionalGroupConfig]] = [],
        manually_excluded_fgs_reactant2: Optional[List[FunctionalGroupConfig]] = [],
        branch_step: Optional[Tuple] = None,
        num_reactants: int = None,
        xviewer_connection: XViewer_Connection = None,
    ):
        self.rxn_name = rxn_name
        self.smarts = smarts
        self.associated_validation_assays_reactant_1 = (
            associated_validation_assays_reactant_1
        )
        self.associated_validation_assays_reactant_2 = (
            associated_validation_assays_reactant_2
        )
        self.example_reactant1 = example_reactant1
        self.example_reactant2 = example_reactant2
        self.example_product = example_product
        self.allowed_fgs_reactant1 = allowed_fgs_reactant1
        self.allowed_fgs_reactant2 = allowed_fgs_reactant2
        self.forbidden_fgs_reactant1 = forbidden_fgs_reactant1
        self.forbidden_fgs_reactant2 = forbidden_fgs_reactant2
        self.product_fgs = product_fgs
        self.manually_added_fgs_reactant2 = manually_added_fgs_reactant2
        self.manually_excluded_fgs_reactant2 = manually_excluded_fgs_reactant2
        self.intermediate_isolated = intermediate_isolated
        self.branch_step = branch_step
        self.xviewer_connection = xviewer_connection
        #Some processing is done to calculate an accurate number of reactants
        self.num_reactants = num_reactants
        self.num_reactants = self.smarts.count(".") + 1
        if self.smarts.count(".") > 0 and (
            self.smarts[0] == "(" or self.smarts[-1] == ")"
        ):
            #SMARTS contains component level grouping
            self.num_reactants -= 1
        if "XLIB0012" in self.rxn_name:
            #Workaround, because there is a "." character in one reactant of the reaction string
            self.num_reactants = 2
        if self.num_reactants > 1 and not self.intermediate_isolated:
            #Add everything that was an allowed functional group in the first reactant as a forbidden functional group in the second reactant. Otherwise, the product can react again!
            for fg in self.allowed_fgs_reactant1:
                if fg.fg_name not in [
                    fg2.fg_name for fg2 in self.forbidden_fgs_reactant2
                ] and fg.fg_name not in [
                    fg2.fg_name for fg2 in self.allowed_fgs_reactant2
                ]:
                    self.forbidden_fgs_reactant2.append(fg)


@dataclass
class ValidationStepConfig(Config):
    """
    Attributes:
    ___________

    validation_step_idx:
        -int, the index of the validation step in the parent validation sequence.
    assay_names:
        -List[str], the list of validation assay_names to be looked up
    validation_cutoff:
        -float, the validation cutoff to be applied. Any BBs which have a validation yield less than this would be excluded.
    assay_to_use_when_equal:
        -str, which of the two validation assays to use when the validation yields of assay1 and assay2 are equal
    include_unvalidated_bbs:
        -bool, whether to consider BBs without validation data for inclusion in the BB list
    xviewer_connection:
        -XViewer_Connection, the connection to the XViewer database
    """

    def __init__(
        self,
        validation_step_idx: int,
        assay_names: List[str],
        validation_cutoff: float,
        assay_to_use_when_equal: Optional[str] = None,
        include_unvalidated_bbs: bool = False,
        xviewer_connection: XViewer_Connection = None,
    ):
        self.validation_step_idx = validation_step_idx
        self.assay_names = assay_names
        self.validation_cutoff = validation_cutoff
        self.assay_to_use_when_equal = assay_to_use_when_equal
        self.include_unvalidated_bbs = include_unvalidated_bbs
        self.xviewer_connection = xviewer_connection


@dataclass
class ReactionBranchConfig(Config):
    """
    Attributes:
    ___________
    reaction_branch_idx:
        -int, the index of the reaction branch in the parent cycle
    reaction_sequence:
        -Tuple[ReactionConfig], the (ordered) sequence of reactions to be applied
    validation_sequence:
        -Tuple[ValidationStepConfig], the (ordered) sequence of validation assays to be applied for this reaction branch. Each element in the tuple represents a unique validation sequence
    monosynthon_exemplar:
        -Optional[str], the SMILES string of an example placeholder structure that is capable of reacting in this reaction branch
    xviewer_connection:
        -XViewer_Connection, the connection to the XViewer database
    """

    def __init__(
        self,
        reaction_branch_idx: int,
        reaction_sequence: Tuple[ReactionConfig],
        validation_sequence: Tuple[ValidationStepConfig] = None,
        monosynthon_exemplar: Optional[str] = None,
        xviewer_connection: XViewer_Connection = None,
    ):
        self.reaction_branch_idx = reaction_branch_idx
        self.reaction_sequence = reaction_sequence
        self.validation_sequence = validation_sequence
        self.monosynthon_exemplar = monosynthon_exemplar
        self.xviewer_connection = xviewer_connection


from configs.xsmiles import (
    ARTIFICIAL_ISOTOPES,
    CONNECTIVITY_HANDLES,
    NULL_PATTERN,
)
from configs.stereochemistry import MAX_N_STEREOISOMERS_ALLOWED
from configs.property_profiles import MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5
from enumeration.exhaustive_enumeration import parallel_map_concatenate_xsmiles
from enumeration.monosynthon_enumeration import (
    monosynthon_enum_cycle,
    monosynthon_react_cycle,
)
from filters.classnames import filter_undesirable_classnames
from filters.deck import filter_deck
from filters.duplicate_enums import filter_duplicate_monosynthon_enumerations
from filters.invalid_xsmiles import filter_df_for_invalid_xsmiles
from filters.substructures import filter_reactivity, filter_undesirable_substructures
from filters.tbbs import filter_tbbs
from filters.validation import (
    filter_validation_for_virtual_libraries,
)
from filters.duplicate_xbbids import filter_duplicate_xbbids
from filters.stereoisomers import filter_stereoisomers
from filters.mw import filter_mw
from filters.virtual_bbs import filter_virtual_bbs
from tools.count_stereoisomers import parallel_map_calculate_num_stereoisomers
from tools.generate_combinations import (
    generate_unique_combinations_with_at_least_one_instruction,
    generate_combinations_with_no_instruction,
)
from transforms.rdkit_mols import (
    parallel_map_compute_mol_from_smiles,
    compute_mol_from_smiles,
    parallel_map_compute_smiles_from_mol,
    compute_smiles_from_mol,
    parallel_map_canonicalize_smiles_in_chunks
)
from transforms.react import react_bimolecular, react_unimolecular
from transforms.monosynthon_exemplar_builder import build_monosynthon_exemplar
from transforms.xsmiles import (
    set_isotopes,
    reset_isotopes,
    fragment_into_xsmiles_precursor,
    add_connectivity_handles,
    parallel_map_set_isotopes,
    parallel_map_reset_isotopes,
    parallel_map_fragment_into_xsmiles_precursors,
    parallel_map_add_connectivity_handles,
    parallel_map_compute_xsmiles,
)
from xviewer_queries.building_blocks import query_candidate_bbs
from xviewer_queries.functional_groups import (
    query_functional_groups_by_first_reaction_in_reaction_branch,
)
from xviewer_queries.xsmiles import query_xsmiles, query_xsmiles_placeholders
from profiling.calculate_physchem import calculate_properties
from profiling.plot_physchem import plot_property_profiles


@dataclass
class LinkerCycleConfig(Config):
    """
    Attributes:
    ___________

    encoding_tag_set:
        -Union[str, None], the tag set used to encode the linker cycle.
    xbbid_smiles_dict:
        -Dict[str, str] of the form {XBBID:SMILES}
    default_smiles:
        -str, the default SMILES strings to be used in certain enumerations (disynthon, instance)
    df:
        -pd.DataFrame, dataframe containing information about the linker BBs
    monosynthons:
        -List[Chem.rdchem.Mol], the monosynthons for the linker cycle
    bb_xsmiles_precursors:
        -List[str], The precursors to the BB XSMILES. These are valid SMILES strings.
    bb_xsmiles:
        -List[str], The BB XSMILES. These are invalid SMILES strings.
    placeholder_xsmiles_precursor:
        -str, the precursor to the placeholder XSMILES. This is a valid SMILES string.
    placeholder_xsmiles:
        -str, the placeholder XSMILES. This is an invalid SMILES string.
    xviewer_connection:
        -XViewer_Connection, the connection to the XViewer database
    
    Methods:
    ________

    initialize_bb_list():
        Initializes the building block list for the linker cycle by populating self.df
        Args: None
        Returns: None
    
    monosynthon_enumeration()
        Performs monosynthon enumeration on the linker cycle.
        Args:
            cycle_names: Tuple[str]
            placeholder_reaction_branches: Tuple[ReactionBranchConfig]
            placeholder_bb_exemplars: Tuple[str]
            compute_xsmiles: bool
            placeholder_xsmiles: Tuple[str]
        Returns: None

    """

    def __init__(
        self,
        encoding_tag_set: Union[str, None],
        xbbid_smiles_dict: Dict[str, str],
        default_smiles: str,
        df: pd.DataFrame = None,
        monosynthons: List[Chem.rdchem.Mol] = None,
        bb_xsmiles_precursors: List[str] = None,
        bb_xsmiles: List[str] = None,
        placeholder_xsmiles_precursor: str = "",
        placeholder_xsmiles: str = "",
        xviewer_connection: XViewer_Connection = None,
    ):
        if encoding_tag_set in ["", "None"]:
            self.encoding_tag_set = None
        else:
            self.encoding_tag_set = encoding_tag_set
        self.xbbid_smiles_dict = xbbid_smiles_dict
        self.default_smiles = default_smiles
        self.df = df
        self.monosynthons = monosynthons
        self.bb_xsmiles_precursors = bb_xsmiles_precursors
        self.bb_xsmiles = bb_xsmiles
        self.placeholder_xsmiles_precursor = placeholder_xsmiles_precursor
        self.placeholder_xsmiles = placeholder_xsmiles
        self.xviewer_connection = xviewer_connection

    def initialize_bb_list(self) -> None:
        df = query_candidate_bbs(xviewer_connection=self.xviewer_connection)

        # Filter for only BBs used in the linker cycle
        df = df[df["XBBID"].isin(list(self.xbbid_smiles_dict.keys()))]
        df = df.reset_index(drop=True)

        df["Encoding Tag Set"] = [self.encoding_tag_set] * len(df)
        df["LOCKSMILES"] = [self.xbbid_smiles_dict[x] for x in list(df["XBBID"])]
        assert len(df) > 0
        self.df = df
        return None

    def monosynthon_enumeration(
        self,
        cycle_names: Tuple[str],
        placeholder_reaction_branches: Tuple[ReactionBranchConfig],
        placeholder_bb_exemplars: Tuple[str],
        compute_xsmiles: bool = False,
        placeholder_xsmiles: Tuple[str] = None,
    ) -> None:

        assert (
            len(cycle_names)
            == len(placeholder_reaction_branches)
            == len(placeholder_bb_exemplars)
        )

        self.monosynthons = parallel_map_compute_mol_from_smiles(
            list(self.df["LOCKSMILES"])
        )
        self.monosynthons = parallel_map_set_isotopes(
            self.monosynthons, ARTIFICIAL_ISOTOPES["L"]
        )

        for i in range(len(cycle_names)):
            cycle_name = cycle_names[i]
            print(f"Reacting cycle {cycle_name} for cycle L monosynthon enumeration.")
            self.df, self.monosynthons = monosynthon_react_cycle(
                self.df,
                cycle_name,
                placeholder_reaction_branches[i].reaction_sequence,
                compute_mol_from_smiles(placeholder_bb_exemplars[i]),
                self.monosynthons,
            )

            self.monosynthons = parallel_map_set_isotopes(
                self.monosynthons, ARTIFICIAL_ISOTOPES[cycle_name]
            )

        self.df.rename(
            columns={"SMILES after " + cycle_names[-1]: "Monosynthon SMILES"},
            inplace=True,
        )

        if compute_xsmiles:
            self.bb_xsmiles_precursors = parallel_map_fragment_into_xsmiles_precursors(
                self.monosynthons, "L"
            )
            self.bb_xsmiles = parallel_map_add_connectivity_handles(
                self.bb_xsmiles_precursors
            )
            self.bb_xsmiles = [x + "." for x in self.bb_xsmiles]

            # Replace any null patterns, if it is unambiguous. If ambiguous, the null patterns will remain present.
            null_replacement_pattern = CONNECTIVITY_HANDLES[NULL_PATTERN]
            for xsmiles in self.bb_xsmiles:
                if (
                    xsmiles.count("%") == 1
                    and xsmiles.count(null_replacement_pattern) == 0
                ):  # A non-blank building block, with exactly one connection
                    null_replacement_pattern = (
                        "%" + xsmiles.split("%")[1][:2]
                    )  # use the same connectivity handle that was present in the non-blank building block
                    break
            self.bb_xsmiles = [
                self.bb_xsmiles[i].replace(
                    CONNECTIVITY_HANDLES[NULL_PATTERN], null_replacement_pattern
                )
                for i in range(len(self.bb_xsmiles))
            ]

            # Remove the linker portion from the first cycle XSMILES placeholder string
            placeholder_xsmiles = list(placeholder_xsmiles)
            placeholder_xsmiles[0] = placeholder_xsmiles[0].split(".")[1] + "."
            placeholder_xsmiles = tuple(placeholder_xsmiles)

            # Compute XSMILES
            xsmiles_fragment_lists = [self.bb_xsmiles]  # Will become a list of lists
            xsmiles_combo_list = []  # Will become a list of tuples
            for i in range(len(cycle_names)):
                xsmiles_fragment_lists.append(
                    [placeholder_xsmiles[i]] * len(self.bb_xsmiles)
                )

            for i in range(len(self.bb_xsmiles)):
                xsmiles_fragments = []
                for xsmiles_fragment_list in xsmiles_fragment_lists:
                    xsmiles_fragments.append(xsmiles_fragment_list[i])
                xsmiles_combo_list.append(tuple(xsmiles_fragments))
            monosynthon_xsmiles_list = parallel_map_compute_xsmiles(xsmiles_combo_list)

            self.df["BB XSMILES Precursors"] = self.bb_xsmiles_precursors
            self.df["BB XSMILES"] = self.bb_xsmiles
            self.df["Monosynthon XSMILES"] = monosynthon_xsmiles_list
            canonical_xsmiles_list = parallel_map_compute_smiles_from_mol(
                parallel_map_compute_mol_from_smiles(monosynthon_xsmiles_list)
            )
            self.df["Monosynthon Canonical XSMILES"] = canonical_xsmiles_list
            self.df["RXN-XSMILES Match"] = [
                list(self.df["Monosynthon SMILES"])[i] == canonical_xsmiles_list[i]
                for i in range(len(self.bb_xsmiles))
            ]
            self.df["XSMILES Valid"] = [
                list(self.df["RXN-XSMILES Match"])[i]
                and canonical_xsmiles_list[i] != ""
                for i in range(len(self.bb_xsmiles))
            ]

        self.monosynthons = parallel_map_reset_isotopes(self.monosynthons)
        self.df["Monosynthon MW"] = [MolWt(m) for m in self.monosynthons]

        return None


@dataclass
class CycleConfig(Config):
    """
    Attributes:
    ___________

    cycle_name:
        -str, the name of the cycle (e.g. "B")
    library_name:
        -str, the name of the library (e.g. "XLIB0001")
    reaction_branches:
        -Tuple[ReactionBranchConfig], the collection of reaction branches that make up the synthetic scheme of this cycle
    working_df:
        -pd.DataFrame, the current list of building blocks which are considered to be suitable for library synthesis, based on filters which it has passed through so far.
    filtered_OUT_df:
        -pd.DataFrame, the list of building blocks which failed to pass at least one filter, and are therefore not considered suitable for library synthesis.
    split_size:
        -int, the number of building blocks used in this cycle, defaults to 0
    filter_for_deck:
        -bool, whether to include only BBs assigned to the library's deck, defaults to False
    bbs:
        -List[Chem.rdchem.Mol], the RDKit mol representations of the SMILES from the working_df
    monosynthons:
        -List[Chem.rdchem.Mol], the enumerated monosynthons for this cycle
    bb_xsmiles_precursors:
        -List[str], the precursor to the XSMILES strings for the BBs at this cycle
    bb_xsmiles:
        -List[str], the XSMILES strings for the BBs at this cycle
    forced_inclusion_bb_dict:
        -Dict, of the form {XBBID: reaction_branch_name}. Only building blocks with XBBIDs on this list will be considered for subsequent processing, and will be used in the specified reaction branch
    forced_exclusion_bbs:
        -List[str], building blocks with these XBBIDs will be manually filtered out of the list (for every reaction branch)
    placeholder_xsmiles_precursor:
        -str, the precursor to the placeholder_xsmiles, with connectivity handles replaced by dummy atoms (so, does correspond to a valid molecule)
    placeholder_xsmiles:
        -str, the XSMILES to use as a placeholder for this cycle, when enumerating (does not correspond to a valid molecule on its own)
    perform_diversity_selection:
        -bool, whether to perform a diversity selection to filter down to a target number of BBs
    target_number_of_bbs:
        -int, the target number of BBs to filter down to using a diversity selection
    xsmiles_dict_no_instruction:
        -dict of the form {XBBID (str): 2-Tuple[XSMILES (str), INSTRUCTION (str)]}, for use in XSMILES disynthon or instance enumeration
    xsmiles_dict_with_instruction:
        -dict of the form {XBBID (str): 2-Tuple[XSMILES (str), INSTRUCTION (str)]}, for use in XSMILES disynthon or instance enumeration
    xviewer_connection:
        -XViewer_Connection, the connection to the XViewer database

    Methods:
    ________

    update_working_df():
        Moves building blocks which did not pass at least one filter from self.working_df to self.filtered_OUT_df and updates the cycle split size.
        Args: None
        Returns: None
    
    initialize_bb_list():
        Populates self.working_df with candidate BBs for this cycle
        Args: None
        Returns: None

    query_compatible_bbs():
        Filters self.working_df for only BBs which are compatible with the reaction scheme
        Args: None
        Returns: None
    
    convert_smiles_to_mols():
        Populates self.bbs by converting the SMILES column in self.working_df to RDKit mols
        Args: None
        Returns: None

    filter_classnames():
        Filters self.working_df for only BBs without undesirable classnames
        Args: None
        Returns: None

    filter_duplicate_enums():
        Filters self.working_df for only BBs with unique monosynthon enumeration structures
        Args: None
        Returns: None

    filter_reactivity():
        Filters self.working_df for only BBs which have fundamentally compatible reactivity for their reaction branch
        Args:
            -override_multiple_substructure_matches, bool
        Returns: None
     
    filter_undesirable_substructures():
        Filters self.working_df for only BBs which do not have undesirable substructures. "Alerts" are filtered out, "Warns" are flagged.
        Args:
            -fgs_to_search, List[FunctionalGroupConfig]
        Returns: None

    filter_tbbs():
        Filters self.working_df for BBs without a "TBB" prefix in their IDs
        Args: None
        Returns: None

    filter_virtual_bbs():
        Filters self.working_df for only BBs which are not virtual
        Args: None
        Returns: None

    filter_validation_for_virtual_libraries():
        Filters self.working_df for only BBs which pass validation. Used only for virtual libraries.
        Args: None
        Returns: None
    
    filter_duplicate_ids():
        Filters self.working_df for only BBs with unique IDs
        Args: None
        Returns: None

    build_monosynthon_exemplar():
        Populates self.bb_exemplar
        Args: None
        Returns: None

    calculate_num_stereoisomers():
        Adds a column corresponding to the number of stereoisomers in the monosynthons for this cycle to self.working_df
        Args:
            -col_name_of_monosynthon_smiles, str
        Returns: None

    monosynthon_enumeration():
        Takes the BBs of this cycle through the entire reaction sequence of the library, using placeholders for all other cycles
        Args:
            -linkersmiles, str
            -cycle_names, Tuple[str]
            -placeholder_reaction_branches, Tuple[ReactionBranchConfig]
            -compute_xsmiles, bool
            -placeholder_xsmiles, Tuple[str]
            -linker_cycle_encoded, bool
            -linker_placeholder_xsmiles, str
        Returns: None

    filter_mw():
        Filters self.working_df for only monosynthons with MW less than or equal to `mw`
        Args:
            -mw, int
        Returns: None
    
    filter_stereoisomers():
        Filters self.working_df for only monosynthons with no more than the number of stereoisomers set as `stereoisomer_threshold`
        Args:
            -stereoisomer_threshold, int
        Returns: None
    
    filter_cycle_for_deck():
        If self.filter_for_deck = True for this cycle, will filter self.working_df for only BBs assigned to the same deck set as `deck`
        Args:
            -deck, str
        Returns: None

    check_bb_count():
        If the user suppied a BB list for this cycle, ensures that the same number of BBs which they supplied is represented in either self.working_df or self.filtered_OUT_df
        Args: None
        Returns: None
    
    lookup_bb_xsmiles():
        Looks up the BB XSMILES for this cycle
        Args:
            -xsmiles_file, str
            -use_int_bbids, bool
        Returns: None
    
    lookup_placeholder_xsmiles():
        Looks up the placeholder XSMILES for this cycle
        Args:
            -xsmiles_placeholder_file, str
        Returns: None        

    filter_invalid_xsmiles():
        Filters self.working_df for only BBs with valid XSMILES. Used primarily for virtual library enumeration
        Args: None
        Returns: None
    
    """

    def __init__(
        self,
        cycle_name: str,
        library_name: str,
        reaction_branches: Tuple[ReactionBranchConfig],
        working_df: pd.DataFrame = None,
        filtered_OUT_df: pd.DataFrame = None,
        split_size: int = None,
        filter_for_deck: bool = False,
        bbs: List[Chem.rdchem.Mol] = None,
        monosynthons: List[Chem.rdchem.Mol] = None,
        bb_xsmiles_precursors: List[str] = None,
        bb_xsmiles: List[str] = None,
        forced_inclusion_bb_dict: Dict = None,
        forced_exclusion_bbs: List[str] = None,
        starting_materials: List[Chem.rdchem.Mol] = None,
        products: List[Chem.rdchem.Mol] = None,
        placeholder_xsmiles_precursor: str = None,
        placeholder_xsmiles: str = None,
        perform_diversity_selection: bool = False,
        target_number_of_bbs: int = None,
        xsmiles_dict_no_instruction: dict = None,
        xsmiles_dict_with_instruction: dict = None,
        xviewer_connection: XViewer_Connection = None,
    ):
        self.cycle_name = cycle_name
        self.library_name = library_name
        self.reaction_branches = reaction_branches
        self.working_df = working_df
        self.filtered_OUT_df = filtered_OUT_df
        self.split_size = split_size
        self.filter_for_deck = filter_for_deck
        self.bbs = bbs
        self.monosynthons = monosynthons
        self.bb_xsmiles_precursors = bb_xsmiles_precursors
        self.bb_xsmiles = bb_xsmiles
        self.forced_inclusion_bb_dict = forced_inclusion_bb_dict
        self.forced_exclusion_bbs = forced_exclusion_bbs
        self.starting_materials = starting_materials
        self.products = products
        self.placeholder_xsmiles_precursor = placeholder_xsmiles_precursor
        self.placeholder_xsmiles = placeholder_xsmiles
        self.perform_diversity_selection = perform_diversity_selection
        self.target_number_of_bbs = target_number_of_bbs
        self.xsmiles_dict_no_instruction = xsmiles_dict_no_instruction
        self.xsmiles_dict_with_instruction = xsmiles_dict_with_instruction
        self.xviewer_connection = xviewer_connection

    def update_working_df(self) -> None:
        keep = self.working_df[self.working_df["Include"] == "Yes"]
        filter_out = self.working_df[self.working_df["Include"] != "Yes"]
        if self.bbs is not None:
            self.bbs = [
                self.bbs[i]
                for i in list(self.working_df.index.values)
                if i in list(keep.index.values)
            ]
        if self.monosynthons is not None:
            self.monosynthons = [
                self.monosynthons[i]
                for i in list(self.working_df.index.values)
                if i in list(keep.index.values)
            ]
        self.working_df = keep.reset_index(drop=True)
        self.split_size = len(self.working_df)
        if self.filtered_OUT_df is None:
            self.filtered_OUT_df = filter_out
        else:
            new_columns = [
                col
                for col in list(self.working_df.columns)
                if col not in list(self.filtered_OUT_df.columns)
            ]
            for col in new_columns:
                if col == "XSMILES":
                    self.filtered_OUT_df[col] = [""] * len(self.filtered_OUT_df)
                else:
                    self.filtered_OUT_df[col] = ["Not tested"] * len(
                        self.filtered_OUT_df
                    )
            self.filtered_OUT_df = pd.concat(
                [self.filtered_OUT_df, filter_out]
            ).reset_index(drop=True)

        # If we have no BBs remaining, then write out the filtered out BBs for inspection, and exit.
        if self.split_size == 0:
            filepath = f"{self.library_name}_cyc_{self.cycle_name}_filtered_out_BBs.csv"
            self.filtered_OUT_df.to_csv(filepath, index=False)
            raise Exception(
                f"Cycle {self.cycle_name} has no remaining BBs! Wrote filtered out BBs to ./{filepath}"
            )

        # If one of the non-blank reaction branches has no remaining BBs, print a warning
        currently_present_reaction_branches = list(
            set(list(self.working_df["Reaction Sequence ID"]))
        )
        if len(currently_present_reaction_branches) != len(self.reaction_branches):
            not_represented_reaction_branches = [
                rb.reaction_branch_idx
                for rb in self.reaction_branches
                if rb.reaction_branch_idx not in currently_present_reaction_branches
            ]
            print(
                f"\nWarning: no remaining BBs in cycle {self.cycle_name} reaction branch(es) {not_represented_reaction_branches}!\n"
            )

        return None

    def initialize_bb_list(self) -> None:
        df = query_candidate_bbs(xviewer_connection=self.xviewer_connection)
        max_num_reactions = max(
            [len(branch.reaction_sequence) for branch in self.reaction_branches]
        )

        # Create empty columns which will eventually describe the reaction sequence for each BB
        df["Reaction Sequence ID"] = [""] * len(df)
        for i in range(max_num_reactions):
            df["Reaction " + str(i + 1)] = [""] * len(df)

        # Filter the df to only include BBs which were present in the input list. This also populates the reaction sequence columns
        if self.forced_inclusion_bb_dict is not None:
            reaction_branch_dict = {
                branch.reaction_branch_idx: branch for branch in self.reaction_branches
            }
            mask = []
            for i in range(len(df)):
                mask.append(
                    df["XBBID"][i] in list(self.forced_inclusion_bb_dict.keys())
                )
            df = df[mask]
            df = df.reset_index(drop=True)
            for i in range(len(df)):
                xbbid = df["XBBID"][i]
                branch = reaction_branch_dict[self.forced_inclusion_bb_dict[xbbid]]
                for j in range(len(branch.reaction_sequence)):
                    df.loc[i, ["Reaction Sequence ID"]] = branch.reaction_branch_idx
                    df.loc[i, ["Reaction " + str(j + 1)]] = branch.reaction_sequence[
                        j
                    ].rxn_name

        self.working_df = df
        self.split_size = len(self.working_df)
        self.bbs = parallel_map_compute_mol_from_smiles(list(self.working_df["SMILES"]))
        return None

    def query_compatible_bbs(self, is_postsynthesis: bool = False) -> None:
        # if we had already populated the reaction sequence columns because of forced inclusions, we don't want to overwrite it. Otherwise, we still need to populate them
        overwrite_existing_reaction_assignments = self.forced_inclusion_bb_dict is None
        df_list = [
            query_functional_groups_by_first_reaction_in_reaction_branch(
                reaction_branch,
                self.working_df,
                self.bbs,
                overwrite_existing_reaction_assignments,
                self.forced_inclusion_bb_dict,
            )
            for reaction_branch in self.reaction_branches
        ]
        matching_df = pd.concat(df_list)
        matching_df = matching_df.reset_index(drop=True)
        if len(matching_df) == 0:
            raise Exception(
                f"No matching BBs were found for cycle {self.cycle_name} after the SMARTS queries were applied"
            )
        if self.forced_inclusion_bb_dict is not None:
            unmatching_xbbids = []
            for i in range(len(self.working_df)):
                xbbid = self.working_df["XBBID"][i]
                if xbbid not in list(matching_df["XBBID"]):
                    unmatching_xbbids.append(xbbid)
            if len(unmatching_xbbids) > 0:
                mask = []
                for i in range(len(self.working_df)):
                    mask.append(self.working_df["XBBID"][i] in unmatching_xbbids)
                self.filtered_OUT_df = self.working_df[mask]
                self.filtered_OUT_df = self.filtered_OUT_df.reset_index(drop=True)
                self.filtered_OUT_df["Include"] = ["No"] * len(self.filtered_OUT_df)
                self.filtered_OUT_df["Contains All Compatible FGs"] = ["No"] * len(
                    self.filtered_OUT_df
                )
        self.working_df = matching_df.sort_values(
            by="Reaction Sequence ID", ascending=True
        ).reset_index(drop=True)
        if is_postsynthesis:
            # We need to append the filtered out BBs to the working_df, so that we can filter them out later
            self.working_df = (
                pd.concat([self.working_df, self.filtered_OUT_df])
                .sort_values(by="Reaction Sequence ID", ascending=True)
                .reset_index(drop=True)
            )
            self.filtered_OUT_df = None
            self.bbs = parallel_map_compute_mol_from_smiles(
                list(self.working_df["SMILES"])
            )
        else:
            self.bbs = parallel_map_compute_mol_from_smiles(
                list(self.working_df["SMILES"])
            )
            self.update_working_df()

        # Now exclude any BBs which the user wanted to force exclusion for
        if self.forced_exclusion_bbs is not None:
            df = self.working_df
            for i in range(len(df)):
                if df["XBBID"][i] in self.forced_exclusion_bbs:
                    df.loc[i, "Excluded By User (XBBID)"] = "Yes"
                    df.loc[i, "Include"] = "No"
            if not is_postsynthesis:
                self.update_working_df()
        return None

    def convert_smiles_to_mols(self) -> None:
        self.bbs = parallel_map_compute_mol_from_smiles(list(self.working_df["SMILES"]))
        return None

    def filter_classnames(self) -> None:
        self.working_df = filter_undesirable_classnames(self.working_df)
        self.update_working_df()
        return None

    def filter_duplicate_enums(self) -> None:
        self.working_df = filter_duplicate_monosynthon_enumerations(
            self.working_df, self.monosynthons
        )
        self.update_working_df()
        return None

    def filter_reactivity(
        self, override_multiple_substructure_matches: bool = False
    ) -> None:
        df_chunks = []
        for branch in self.reaction_branches:
            branch_df = self.working_df[
                self.working_df["Reaction Sequence ID"] == branch.reaction_branch_idx
            ]
            bbs = [
                self.bbs[i]
                for i in range(len(self.bbs))
                if i in list(branch_df.index.values)
            ]
            branch_df = branch_df.reset_index(drop=True)
            df_chunks.append(
                filter_reactivity(
                    branch_df,
                    bbs,
                    branch.reaction_sequence,
                    override_multiple_substructure_matches=override_multiple_substructure_matches,
                )
            )
        self.working_df = pd.concat(df_chunks).reset_index(drop=True)
        self.bbs = parallel_map_compute_mol_from_smiles(list(self.working_df["SMILES"]))
        self.update_working_df()
        return None

    def filter_undesirable_substructures(self, fgs_to_search: List[FunctionalGroupConfig]) -> None:
        self.working_df = filter_undesirable_substructures(
            self.working_df, self.monosynthons, fgs_to_search
        )
        self.update_working_df()
        return None

    def filter_tbbs(self) -> None:
        self.working_df = filter_tbbs(self.working_df)
        self.update_working_df()
        return None

    def filter_virtual_bbs(self) -> None:
        self.working_df = filter_virtual_bbs(self.working_df)
        self.update_working_df()
        return None

    def filter_validation_for_virtual_libraries(self) -> None:
        self.working_df = filter_validation_for_virtual_libraries(self.working_df)
        self.update_working_df()
        return None

    def filter_duplicate_ids(self) -> None:
        self.working_df = filter_duplicate_xbbids(self.working_df)
        self.update_working_df()
        return None

    def build_monosynthon_exemplar(self) -> None:
        self.bb_exemplar = build_monosynthon_exemplar(self.reaction_branches)
        return None

    def calculate_num_stereoisomers(self, col_name_of_monosynthon_smiles: str) -> None:
        df = self.working_df

        # RDKit's EnumerateStereoisomers method does not work on molecule lists generated with parallel map for some reason, but it does work on molecules that were generated
        # individually, hence the workaround below where we regenerate the monosynthons one by one.
        monosynthons = [
            compute_mol_from_smiles(smi)
            for smi in list(df[col_name_of_monosynthon_smiles])
        ]

        df["Number of Stereoisomers in Monosynthon"] = (
            parallel_map_calculate_num_stereoisomers(
                monosynthons, list(df["STEREOCHEMISTRYNAME"])
            )
        )
        self.working_df = df
        return None

    def monosynthon_enumeration(
        self,
        linkersmiles: str,
        cycle_names: Tuple[str],
        placeholder_reaction_branches: Tuple[ReactionBranchConfig],
        placeholder_bb_exemplars: Tuple[str],
        compute_xsmiles: bool = False,
        placeholder_xsmiles: Tuple[str] = None,
        linker_cycle_encoded: bool = False,
        linker_placeholder_xsmiles: str = None,
    ) -> None:

        assert (
            len(cycle_names)
            == len(placeholder_reaction_branches)
            == len(placeholder_bb_exemplars)
        )
        self.monosynthons = [compute_mol_from_smiles(linkersmiles)] * self.split_size
        self.monosynthons = parallel_map_set_isotopes(
            self.monosynthons, ARTIFICIAL_ISOTOPES["L"]
        )

        for i in range(len(cycle_names)):
            cycle_name = cycle_names[i]
            if self.cycle_name == cycle_name:
                df_chunks = []
                monosynthon_chunks = []
                bb_chunks = []
                for branch in self.reaction_branches:
                    branch_df = self.working_df[
                        self.working_df["Reaction Sequence ID"]
                        == branch.reaction_branch_idx
                    ]
                    if len(branch_df) > 0:
                        print(
                            f"Enumerating cycle {cycle_name} ({self.split_size} instances) reaction branch {branch.reaction_branch_idx} ({len(branch_df)} instances) for monosynthon enumeration"
                        )
                        bbs = [
                            self.bbs[x]
                            for x in range(len(self.bbs))
                            if x in list(branch_df.index.values)
                        ]
                        monosynthons = [
                            self.monosynthons[x]
                            for x in range(len(self.monosynthons))
                            if x in list(branch_df.index.values)
                        ]
                        branch_df = branch_df.reset_index(drop=True)
                        new_df, new_monosynthons = monosynthon_enum_cycle(
                            self.cycle_name,
                            branch.reaction_sequence,
                            branch_df,
                            bbs,
                            monosynthons,
                        )
                        df_chunks.append(new_df)
                        monosynthon_chunks.extend(new_monosynthons)
                        bb_chunks.extend(bbs)
                self.working_df = pd.concat(df_chunks).reset_index(drop=True)
                self.monosynthons = monosynthon_chunks
                self.bbs = bb_chunks

            else:
                print(
                    f"Reacting cycle {cycle_name} for cycle {self.cycle_name} monosynthon enumeration ({self.split_size} instances)"
                )
                self.working_df, self.monosynthons = monosynthon_react_cycle(
                    self.working_df,
                    cycle_name,
                    placeholder_reaction_branches[i].reaction_sequence,
                    compute_mol_from_smiles(placeholder_bb_exemplars[i]),
                    self.monosynthons,
                )

            self.monosynthons = parallel_map_set_isotopes(
                self.monosynthons, ARTIFICIAL_ISOTOPES[cycle_name]
            )

        self.working_df.rename(
            columns={"SMILES after " + cycle_names[-1]: "Monosynthon SMILES"},
            inplace=True,
        )

        if compute_xsmiles:
            linker_xsmiles_precursors = parallel_map_fragment_into_xsmiles_precursors(
                self.monosynthons, "L"
            )
            linker_xsmiles = parallel_map_add_connectivity_handles(
                linker_xsmiles_precursors
            )

            self.bb_xsmiles_precursors = parallel_map_fragment_into_xsmiles_precursors(
                self.monosynthons, self.cycle_name
            )
            self.bb_xsmiles = parallel_map_add_connectivity_handles(
                self.bb_xsmiles_precursors
            )

            # Replace any null patterns, if it is unambiguous. If ambiguous, the null patterns will remain present.
            null_replacement_pattern = CONNECTIVITY_HANDLES[NULL_PATTERN]
            for xsmiles in self.bb_xsmiles:
                if (
                    xsmiles.count("%") == 1
                    and xsmiles.count(null_replacement_pattern) == 0
                ):  # A non-blank building block, with exactly one connection
                    null_replacement_pattern = (
                        "%" + xsmiles.split("%")[1][:2]
                    )  # use the same connectivity handle that was present in the non-blank building block
                    break
            self.bb_xsmiles = [
                self.bb_xsmiles[i].replace(
                    CONNECTIVITY_HANDLES[NULL_PATTERN], null_replacement_pattern
                )
                for i in range(len(self.bb_xsmiles))
            ]

            if linker_cycle_encoded:
                xsmiles_fragment_lists = [
                    [linker_placeholder_xsmiles] * len(self.bb_xsmiles)
                ]  # Will become a list of lists
            else:
                # If this is the first cycle, generate and concatenate the linker XSMILES to this cycle's XSMILES
                if self.cycle_name == cycle_names[0]:
                    self.bb_xsmiles = [
                        linker_xsmiles[i] + "." + self.bb_xsmiles[i]
                        for i in range(len(self.bb_xsmiles))
                    ]
                    # The enum now needs to not start with the linker placeholder SMILES, since we just concatenated it to the BB XSMILES
                    xsmiles_fragment_lists = []  # Will become a list of lists
                else:
                    xsmiles_fragment_lists = [
                        [linker_placeholder_xsmiles] * len(self.bb_xsmiles)
                    ]  # Will become a list of lists

            # Add '.' to the end of the XSMILES as long as it's not the last cycle, so they will be different fragments upon string concatenation with other cycles
            if self.cycle_name != cycle_names[-1]:
                self.bb_xsmiles = [s + "." for s in self.bb_xsmiles]

            # Remove the linker portion from the first cycle XSMILES placeholder string
            placeholder_xsmiles = list(placeholder_xsmiles)
            placeholder_xsmiles[0] = placeholder_xsmiles[0].split(".")[1] + "."
            placeholder_xsmiles = tuple(placeholder_xsmiles)

            # Compute XSMILES
            xsmiles_combo_list = []  # Will become a list of tuples
            for i in range(len(cycle_names)):
                cycle_name = cycle_names[i]
                if cycle_name != self.cycle_name:
                    xsmiles_fragment_lists.append(
                        [placeholder_xsmiles[i]] * len(self.bb_xsmiles)
                    )
                else:
                    xsmiles_fragment_lists.append(self.bb_xsmiles)
            for i in range(len(self.bb_xsmiles)):
                xsmiles_fragments = []
                for xsmiles_fragment_list in xsmiles_fragment_lists:
                    xsmiles_fragments.append(xsmiles_fragment_list[i])
                xsmiles_combo_list.append(tuple(xsmiles_fragments))
            monosynthon_xsmiles_list = parallel_map_compute_xsmiles(xsmiles_combo_list)

            self.working_df["BB XSMILES Precursors"] = self.bb_xsmiles_precursors
            self.working_df["BB XSMILES"] = self.bb_xsmiles
            self.working_df["Monosynthon XSMILES"] = monosynthon_xsmiles_list
            canonical_xsmiles_list = parallel_map_compute_smiles_from_mol(
                parallel_map_compute_mol_from_smiles(monosynthon_xsmiles_list)
            )
            self.working_df["Monosynthon Canonical XSMILES"] = canonical_xsmiles_list
            self.working_df["RXN-XSMILES Match"] = [
                list(self.working_df["Monosynthon SMILES"])[i]
                == canonical_xsmiles_list[i]
                for i in range(len(self.bb_xsmiles))
            ]
            self.working_df["XSMILES Valid"] = [
                list(self.working_df["RXN-XSMILES Match"])[i]
                and canonical_xsmiles_list[i] != ""
                for i in range(len(self.bb_xsmiles))
            ]

        self.monosynthons = parallel_map_reset_isotopes(self.monosynthons)
        self.working_df["Monosynthon MW"] = [MolWt(m) for m in self.monosynthons]
        self.update_working_df()
        return None

    def filter_mw(self, mw: int = MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5) -> None:
        self.working_df = filter_mw(self.working_df, mw=mw)
        self.update_working_df()
        return None

    def filter_stereoisomers(self, stereoisomer_threshold: int = MAX_N_STEREOISOMERS_ALLOWED) -> None:
        self.working_df = filter_stereoisomers(self.working_df, stereoisomer_threshold = stereoisomer_threshold)
        self.update_working_df()
        return None

    def filter_cycle_for_deck(self, deck: str) -> None:
        if self.filter_for_deck:
            self.working_df = filter_deck(self.working_df, deck)
            self.update_working_df()
        return None

    def check_bb_count(self) -> None:
        if self.forced_inclusion_bb_dict is not None:
            if self.forced_exclusion_bbs is not None:
                expected_bbs = [
                    x
                    for x in list(self.forced_inclusion_bb_dict.keys())
                    if x not in self.forced_exclusion_bbs
                ]
            else:
                expected_bbs = [x for x in list(self.forced_inclusion_bb_dict.keys())]
            found_bbs = list(self.working_df["XBBID"]) + list(
                self.filtered_OUT_df["XBBID"]
            )
            missing_bbs = [x for x in expected_bbs if x not in found_bbs]
            if len(missing_bbs) > 0:
                raise Exception(
                    f"Expected to find the following BBs, but they are missing from the output file: {missing_bbs}"
                )
        return None

    def lookup_bb_xsmiles(
        self, use_int_bbids: bool = False, prefix: str = ''
    ) -> None:
        unique_bb_list_from_input_file = list(set(list(self.working_df["XBBID"])))
        xsmiles_df = query_xsmiles(library_name = self.library_name, cycle_name = self.cycle_name, xviewer_connection = self.xviewer_connection)
        xsmiles_df = xsmiles_df[xsmiles_df["SET_TYPE"] == "DEFAULT"]
        invalid_df = xsmiles_df[xsmiles_df["XSMILES_STATUS"] != "VALID"]
        if len(invalid_df) > 0:
            invalid_bbs = list(invalid_df["XBBID"])
            raise Exception(f"These BBs have invalid XSMILES for cycle {self.cycle_name}: {invalid_bbs}")
        if self.library_name in ['XCLB0130C', 'XCLB0131C'] and self.cycle_name == "D":
            #Add the missing blank BB
            xsmiles_df = pd.DataFrame(
                {
                    'ID':[None],
                    'LIBRARY_NAME':[self.library_name],
                    'CYCLE_NAME':['D'],
                    'SET_ID':[None],
                    'SET_NAME':[f'{self.library_name}_D DEFAULT Set'],
                    'SET_TYPE':['DEFAULT'],
                    'BBID':[860608],
                    'XBBID':['XBB014116'],
                    'XSMILES_ORIG':[''],
                    'XSMILES_FOR_ENUM':[''],
                    'INSTRUCTION':[''],
                    'XSMILES_STATUS':['VALID'],
                }
            )
        xsmiles_df = xsmiles_df.reset_index(drop=True)

        if use_int_bbids:
            xbbid_bbid_dict = {
                xsmiles_df["XBBID"][i]: xsmiles_df["BBID"][i]
                for i in range(len(xsmiles_df))
            }
        else:
            xbbid_bbid_dict = {
                xsmiles_df["XBBID"][i]: xsmiles_df["XBBID"][i]
                for i in range(len(xsmiles_df))
            }

        unique_bb_list_from_xsmiles_file = list(set(list(xsmiles_df["XBBID"])))
        unique_bb_list_from_input_file.sort()
        unique_bb_list_from_xsmiles_file.sort()
        # Make sure no duplicates were present
        try:
            assert len(unique_bb_list_from_input_file) == len(self.working_df)
        except:
            raise Exception(
                f"Length of the UNIQUE BB list ({len(unique_bb_list_from_input_file)}) is different from the length of the input file ({len(self.working_df)})."
            )
        try:
            assert len(unique_bb_list_from_xsmiles_file) == len(xsmiles_df)
        except:
            raise Exception(
                f"Length of the UNIQUE XSMILES records ({len(unique_bb_list_from_xsmiles_file)}) is different from the length of the XSMILES file ({len(xsmiles_df)})."
            )
        # Make sure the BB lists are the exact same
        try:
            assert unique_bb_list_from_input_file == unique_bb_list_from_xsmiles_file
        except:
            raise Exception(
                f"BB lists differ between the input file ({len(unique_bb_list_from_input_file)} records) and the XSMILES file ({len(unique_bb_list_from_xsmiles_file)} records)"
            )

        xsmiles_df_without_instruction = xsmiles_df[xsmiles_df["INSTRUCTION"] == ""]
        xsmiles_df_without_instruction.sort_values(by=["XBBID"], inplace=True)

        xsmiles_df_with_instruction = xsmiles_df[xsmiles_df["INSTRUCTION"] != ""]
        xsmiles_df_with_instruction = xsmiles_df_with_instruction.reset_index(drop=True)
        xsmiles_df_with_instruction = xsmiles_df_with_instruction[
            xsmiles_df_with_instruction["INSTRUCTION"] != "FAIL<>"
        ]  # Remove Failures
        xsmiles_df_with_instruction.sort_values(by=["XBBID"], inplace=True)
        xsmiles_df_without_instruction = xsmiles_df_without_instruction.reset_index(
            drop=True
        )
        xsmiles_df_with_instruction = xsmiles_df_with_instruction.reset_index(drop=True)
        # create dictionaries
        self.xsmiles_dict_no_instruction = {
            xbbid_bbid_dict[xsmiles_df_without_instruction["XBBID"][i]]: (
                prefix + xsmiles_df_without_instruction["XSMILES_FOR_ENUM"][i].strip(),
                xsmiles_df_without_instruction["INSTRUCTION"][i].strip(),
            )
            for i in range(len(xsmiles_df_without_instruction))
        }
        self.xsmiles_dict_with_instruction = {
            xbbid_bbid_dict[xsmiles_df_with_instruction["XBBID"][i]]: (
                prefix + xsmiles_df_with_instruction["XSMILES_FOR_ENUM"][i].strip(),
                xsmiles_df_with_instruction["INSTRUCTION"][i].strip(),
            )
            for i in range(len(xsmiles_df_with_instruction))
        }
        self.split_size = len(self.xsmiles_dict_no_instruction) + len(
            self.xsmiles_dict_with_instruction
        )
        return None

    def lookup_placeholder_xsmiles(self) -> None:
        xsmiles_placeholder_df = query_xsmiles_placeholders(library_name = self.library_name, cycle_name = self.cycle_name, xviewer_connection = self.xviewer_connection)
        if self.library_name in ['XCLB0130C', 'XCLB0131C'] and self.cycle_name == 'D':
            #Add the missing placeholder
            xsmiles_placeholder_df = pd.DataFrame({'PLACEHOLDER_XSMILES':['']})
        xsmiles_placeholder_df = xsmiles_placeholder_df.reset_index(drop=True)
        try:
            assert len(xsmiles_placeholder_df) == 1
        except:
            raise Exception(
                f"The placeholder file for library {self.library_name} cycle {self.cycle_name} has {len(xsmiles_placeholder_df)} rows, it should have 1"
            )
        self.placeholder_xsmiles = str(
            xsmiles_placeholder_df["PLACEHOLDER_XSMILES"][0]
        ).strip()
        return None

    def filter_invalid_xsmiles(self) -> None:
        self.working_df = filter_df_for_invalid_xsmiles(self.working_df)
        self.update_working_df()
        return None


from enumeration.enum_library_subset import enumerate_library_instances


@dataclass
class LibraryConfig(Config):
    """
    Attributes:
    ___________

    lib_name:
        str, the name of the library (e.g. "XLIB0001")
    synonym:
        str, a synonym for the library (e.g. "ODLIB001")
    linker_cycle:
        LinkerCycleConfig, the linker cycle of the library
    cycles:
        Tuple[CycleConfig], the chemistry cycles of the library
    deck (optional):
        str, the deck to which this library belongs (e.g. "DELcore" or "DELflex")
    instance_df:
        pd.DataFrame: contains example instances from the library
    qc_set_df:
        pd.DataFrame: contains the enumeration qc set instances from the library
    instances:
        List[Chem.rdchem.mol]: the instance molecules for the library
    xviewer_connection:
        XViewer_Connection: the connection to the XViewer database

    Methods:
    ________

    enumerate_instances():
        Enumerates instance exemplars for the library (`num_instances` of them)
        Args:
            -num_instances, int
        Returns: None

    enumerate_virtual_instances():
        Enumerates instance exemplars for the virtual library (`num_instances` of them)
        Args:
            -num_instances, int
        Returns: None
    
    enumerate_qc_set():
        Enumerates enough instances to guarantee that every BB is used in at least one instance
        Args: None
        Returns: None
    
    calculate_physchem_properties():
        Calculates physicochemical properties of the instance exemplars for the library
        Args: None
        Returns: None
    
    plot_physchem_property_profiles():
        Exports plots of the physicochemical property profiles
        Args:
            -physchem_property_directory, str
            -keep_autogenerated_property_columns, bool
        Returns: None

    compute_placeholder_xsmiles():
        Generates placeholder XSMILES for each cycle
        Args: None
        Returns: None
    
    validate_placeholder_xsmiles():
        Validates that concatenation of the placeholder XSMILES strings for all cycles generates a valid molecule
        Args: None
        Returns: None
    
    write_excel_files():
        Writes the excel files containing the BB lists for the library
        Args:
            -input_filepath, str
            -output_filepath, str
            -qc_set_filepath, str
            -readidel_exemplar_filepath, str
            -color_list_tabs, bool
            -is_presynthesis, bool
        Returns: None
    
    write_xsmiles_lists():
        Writes the list of BB XSMILES strings for all cycles of the library
        Args:
            -xsmiles_directory, str
        Returns: None
    
    enumerate_disynthons():
        Enumerates disynthons of the specified types for the library ("BC", "BD", ...)
        Args:
            -disynthon_type, Tuple[str]
        Returns:
            pa.Table
    
    enumerate_and_write_instances():
        Enumerates and writes all library instances in chunks to `output_filepath`
        Args:
            -output_filepath, str
        Returns: None
    """

    def __init__(
        self,
        lib_name: str,
        linker_cycle: LinkerCycleConfig,
        cycles: Tuple[CycleConfig],
        synonym:Optional[str] = None,
        deck: Optional[str] = None,
        instance_df: pd.DataFrame = pd.DataFrame(),
        qc_set_df: pd.DataFrame = pd.DataFrame(),
        instances: Optional[List[Chem.rdchem.Mol]] = None,
        xviewer_connection: XViewer_Connection = None,
    ):
        self.lib_name = lib_name
        self.linker_cycle = linker_cycle
        self.cycles = cycles
        self.synonym = synonym
        self.deck = deck
        self.instance_df = instance_df
        self.qc_set_df = qc_set_df
        self.instances = instances
        self.xviewer_connection = xviewer_connection

    def enumerate_instances(self, num_instances: int, add_substructure_warnings:bool = False) -> None:
        self.instance_df, self.instances = enumerate_library_instances(
            self.lib_name,
            self.cycles,
            self.linker_cycle,
            num_instances,
            qc_set=False,
            is_virtual=False,
            add_substructure_warnings=add_substructure_warnings,
        )
        return None

    def enumerate_virtual_instances(self, num_instances: int, add_substructure_warnings:bool = False) -> None:
        self.instance_df, self.instances = enumerate_library_instances(
            self.lib_name,
            self.cycles,
            self.linker_cycle,
            num_instances,
            qc_set=False,
            is_virtual=True,
            add_substructure_warnings = add_substructure_warnings
        )
        return None

    def enumerate_qc_set(self) -> None:
        self.qc_set_df, _ = enumerate_library_instances(
            self.lib_name,
            self.cycles,
            self.linker_cycle,
            qc_set=True,
            is_virtual=False,
        )
        self.qc_set_df["Library"] = [self.lib_name] * len(self.qc_set_df)
        return None

    def calculate_physchem_properties(self) -> None:
        self.instance_df = calculate_properties(self.instance_df, mols=self.instances)
        return None

    def plot_physchem_property_profiles(
        self,
        physchem_profile_directory: str,
        keep_autogenerated_property_columns: bool = False,
    ) -> None:
        self.instance_df = plot_property_profiles(
            self.instance_df,
            physchem_profile_directory,
            keep_autogenerated_property_columns=keep_autogenerated_property_columns,
        )
        return None

    def compute_placeholder_xsmiles(self) -> None:
        placeholder_bbs = [
            compute_mol_from_smiles(cycle.reaction_branches[0].monosynthon_exemplar)
            for cycle in self.cycles
        ]
        for i in range(len(self.cycles)):
            parent_cycle = self.cycles[i]
            placeholder_monosynthon = compute_mol_from_smiles(
                self.linker_cycle.default_smiles
            )
            placeholder_monosynthon = set_isotopes(
                placeholder_monosynthon, ARTIFICIAL_ISOTOPES["L"]
            )
            for j in range(len(self.cycles)):
                current_cycle = self.cycles[j]
                current_cycle_name = current_cycle.cycle_name
                current_reaction_branch = current_cycle.reaction_branches[0]
                bb = deepcopy(placeholder_bbs[j])
                bb = reset_isotopes(bb)
                for rxn in current_reaction_branch.reaction_sequence:
                    if rxn.num_reactants == 2:
                        placeholder_monosynthon = react_bimolecular(
                            (placeholder_monosynthon, bb), rxn.smarts
                        )
                    elif rxn.num_reactants == 1:
                        placeholder_monosynthon = react_unimolecular(
                            placeholder_monosynthon, rxn.smarts
                        )
                    else:
                        raise Exception(
                            f"{rxn.rxn_name} reaction has an unsupported number of reactants, {rxn.num_reactants}"
                        )
                placeholder_monosynthon = set_isotopes(
                    placeholder_monosynthon, ARTIFICIAL_ISOTOPES[current_cycle_name]
                )
            # Now that the full monosynthon has been enumerated, loop over the cycles again to generate the placeholder XSMILES
            placeholder_monosynthon_copy = deepcopy(placeholder_monosynthon)
            xs_precursor = fragment_into_xsmiles_precursor(
                placeholder_monosynthon, parent_cycle.cycle_name
            )
            xs = add_connectivity_handles(xs_precursor)
            linker_xsmiles_precursor = fragment_into_xsmiles_precursor(
                placeholder_monosynthon_copy, "L"
            )
            linker_xsmiles = add_connectivity_handles(linker_xsmiles_precursor) + "."

            if i == 0:
                xs = linker_xsmiles + xs
            if i != range(len(self.cycles))[-1]:
                xs = xs + "."
            setattr(parent_cycle, "placeholder_xsmiles_precursor", xs_precursor)
            setattr(parent_cycle, "placeholder_xsmiles", xs)

        setattr(
            self.linker_cycle, "placeholder_xsmiles_precursor", linker_xsmiles_precursor
        )
        setattr(self.linker_cycle, "placeholder_xsmiles", linker_xsmiles)

        return None

    def validate_placeholder_xsmiles(self) -> None:
        placeholder_xsmiles = ""

        for cycle in self.cycles:
            placeholder_xsmiles = placeholder_xsmiles + cycle.placeholder_xsmiles
        if compute_smiles_from_mol(compute_mol_from_smiles(placeholder_xsmiles)) == "":
            raise Exception(
                f"Placeholder XSMILES validation failed, with placeholder xsmiles: {placeholder_xsmiles}"
            )
        return None

    def write_excel_files(
        self,
        input_filepath: str,
        output_filepath: str,
        qc_set_filepath: str,
        readidel_exemplar_filepath:str,
        color_list_tabs: bool,
        is_presynthesis: bool = False,
    ) -> None:
        xlsx = pd.ExcelFile(input_filepath)

        sheets = {}
        for sheet_name in xlsx.sheet_names:
            if sheet_name in ["Info", "Instance Enumeration", "Enumeration QC Set"] or (
                sheet_name[0] in ["L"] + [cycle.cycle_name for cycle in self.cycles]
                and not sheet_name.startswith("EXAMPLE")
            ):
                sheets[sheet_name] = xlsx.parse(sheet_name)

        sheets["L_lists"] = tuple([self.linker_cycle.df, False])

        for cycle in self.cycles:
            # Tuple: (working_df for this cycle, and bool: whether or not any BBs were filtered out of the list)
            if not is_presynthesis:
                cols_to_exclude = [
                    "Include",
                    "Excluded By User (XBBID)",
                    "Excluded By User (FG)",
                    "Contains All Compatible FGs",
                    "Deck Alert",
                    "Undesirable Classname",
                    "TBB",
                    "Reactive FG Issue",
                    "Failed Enumeration in Cycle",
                    "Incompatible FG",
                    "Capable of Reacting Multiple Times",
                    "Duplicate Monosynthon Enumeration Of",
                    "Alert FGs",
                    "Exclude Due to Alert Substructure",
                    "XBB Used In Another Reaction Sequence",
                    f"More than {str(MAX_N_STEREOISOMERS_ALLOWED)} Stereoisomers in Monosynthon",
                ]
                cycle.working_df = cycle.working_df[
                    [
                        col
                        for col in list(cycle.working_df.columns)
                        if col not in cols_to_exclude
                    ]
                ]
            cycle.working_df = cycle.working_df[
                [col for col in list(cycle.working_df.columns) if col != "Include"]
            ]
            cycle.filtered_OUT_df = cycle.filtered_OUT_df[
                [col for col in list(cycle.filtered_OUT_df.columns) if col != "Include"]
            ]
            
            #Sort values by Reaction Sequence ID and XBBID
            cycle.working_df.sort_values(by=["Reaction Sequence ID", "XBBID"], inplace=True)
            cycle.filtered_OUT_df.sort_values(by=["Reaction Sequence ID", "XBBID"], inplace=True)
            cycle.working_df = cycle.working_df.reset_index(drop=True)
            cycle.filtered_OUT_df = cycle.filtered_OUT_df.reset_index(drop=True)
            
            cycle.filtered_OUT_df.rename(
                columns={
                    name: "EXCLUDED_" + name
                    for name in list(cycle.filtered_OUT_df.columns)
                },
                inplace=True,
            )

            sheets[cycle.cycle_name + "_lists"] = tuple(
                [
                    pd.concat(
                        [
                            cycle.working_df,
                            pd.DataFrame({"": ""}, index=[0]),
                            cycle.filtered_OUT_df,
                        ],
                        axis=1,
                    ),
                    len(cycle.filtered_OUT_df) > 0,
                ]
            )
        sheets["Instance Enumeration"] = self.instance_df
        sheets["Enumeration QC Set"] = self.qc_set_df

        with pd.ExcelWriter(output_filepath, engine="xlsxwriter") as writer:
            workbook = writer.book
            black = workbook.add_format({"bg_color": "#000000"})
            red = workbook.add_format({"bg_color": "#FFC7CE", "font_color": "#9C0006"})
            green = workbook.add_format(
                {"bg_color": "#C6EFCE", "font_color": "#006100"}
            )
            for sheet_name in list(sheets.keys()):
                if "_list" in sheet_name:
                    df = sheets[sheet_name][0].fillna("")
                    contains_filtered_out_bbs = sheets[sheet_name][1]
                else:
                    df = sheets[sheet_name].fillna("")
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                worksheet = writer.sheets[sheet_name]
                coloring_included_xbbid_column = True
                for idx, col in enumerate(df):
                    series = df.iloc[:, idx]
                    width = len(str(series.name)) + 3  # add some extra space
                    worksheet.set_column(idx, idx, width)
                    if "_list" in sheet_name:
                        if str(series.name).endswith("XBBID"):
                            width = max(width, 12)  # Expand XBBID column to fit the length of XBBIDs (9) and some extra space (3)
                            if coloring_included_xbbid_column:
                                worksheet.set_column(idx, idx, width, green)
                                coloring_included_xbbid_column = False
                            else:
                                worksheet.set_column(idx, idx, width, red)
                            if color_list_tabs:
                                if contains_filtered_out_bbs:
                                    worksheet.set_tab_color(
                                        "#9C0006"
                                    )  # Change tab color to red
                                    if "XSMILES Valid" in list(df.columns):
                                        if not all(list(df["XSMILES Valid"])):
                                            worksheet.set_tab_color(
                                                "#7030A0"
                                            )  # Change tab color to purple
                                else:
                                    worksheet.set_tab_color(
                                        "#006100"
                                    )  # Change tab color to green
                                    if "XSMILES Valid" in list(df.columns):
                                        if not all(list(df["XSMILES Valid"])):
                                            worksheet.set_tab_color(
                                                "#00B0F0"
                                            )  # Change tab color to blue
                    if str(series.name).startswith("Unnamed") or str(series.name) == "":
                        worksheet.set_column(idx, idx, 0.3, black)
                if (
                    sheet_name in ["Instance Enumeration", "Enumeration QC Set"]
                    and len(df) > 0
                ):
                    if (
                        list(df["Instance SMILES"]).count("") > 0
                    ):  # If any instances failed to enumerate
                        worksheet.set_tab_color("#9C0006")  # Change tab color to red
                        if "XSMILES Behavior Expected" in list(df.columns):
                            if all(list(df["XSMILES Behavior Expected"])):
                                worksheet.set_tab_color(
                                    "#00B0F0"
                                )  # Change tab color to blue
                    else:
                        worksheet.set_tab_color("#006100")  # Change tab color to green
                        if "XSMILES Behavior Expected" in list(df.columns):
                            if not all(list(df["XSMILES Behavior Expected"])):
                                worksheet.set_tab_color(
                                    "#7030A0"
                                )  # Change tab color to purple
        if len(self.qc_set_df) > 0:
            self.qc_set_df.to_csv(qc_set_filepath, index=False)
        if readidel_exemplar_filepath is not None:
            print(readidel_exemplar_filepath)
            if self.synonym is None:
                raise Exception(f"The library has no synonym, but a readidel exemplar filepath was provided.")
            if not all([f'BB{cycle.cycle_name} EBBID' in list(self.instance_df.columns) for cycle in self.cycles]):
                raise Exception(f'A readidel exemplar filepath was provided, but the instance df does not contain the required columns: {[f"BB{cycle.cycle_name} EBBID" for cycle in self.cycles]}.')
            if len(self.instance_df) < 3000:
                raise Exception(f'The instance df has {len(self.instance_df)} rows, but it must have at least 3000.')
            if any([x == '' for x in list(self.instance_df['Instance SMILES'])]):
                raise Exception(f'The instance df contains empty SMILES strings.')
            for cycle in self.cycles:
                if not all([x.startswith('EBB') for x in list(self.instance_df[f'BB{cycle.cycle_name} EBBID'])]):
                    xbbid_failures = []
                    ebbid_failures = []
                    for i in range(len(self.instance_df)):
                        if not self.instance_df[f'BB{cycle.cycle_name} EBBID'][i].startswith('EBB'):
                            xbbid_failures.append(self.instance_df[f'BB{cycle.cycle_name} ID'][i])
                            ebbid_failures.append(self.instance_df[f'BB{cycle.cycle_name} EBBID'][i])
                    raise Exception(f'The instance df contains EBBIDs which do not start with "EBB". XBBIDs: {xbbid_failures}, EBBIDs: {ebbid_failures}')
            
            readidel_exemplar_df = self.instance_df.head(3000)
            readidel_exemplar_df['Library'] = [self.synonym] * len(readidel_exemplar_df)
            assert len(readidel_exemplar_df) == 3000
            readidel_exemplar_df = readidel_exemplar_df[['Library'] + [f'BB{cycle.cycle_name} EBBID' for cycle in self.cycles] + ['Instance SMILES']]
            column_name_map = {f'BB{cycle.cycle_name} EBBID': f'BB{cycle.cycle_name} ID' for cycle in self.cycles}
            column_name_map.update({'Instance SMILES': 'SMILES'})
            readidel_exemplar_df.rename(columns=column_name_map, inplace=True)
            with pd.ExcelWriter(readidel_exemplar_filepath, engine="xlsxwriter") as writer:
                pandas.io.formats.excel.ExcelFormatter.header_style = None
                readidel_exemplar_df.to_excel(writer, sheet_name="Sheet1", index=False)  
        return None

    def write_xsmiles_lists(self, xsmiles_directory: str) -> None:
        bool_to_string_map = {True: "Y", False: "N", "Not tested": "N"}
        # Handle linker cycle, if it is encoded
        if self.linker_cycle.encoding_tag_set is not None:
            if self.cycles[0].cycle_name == "A":
                raise Exception("Linker cycle is encoded, and there is a cycle A, so the naming of the XSMILES files cannot be determined.")
            data = []

            working_df = self.linker_cycle.df.reset_index(drop=True)
            included_xbbids = list(working_df["XBBID"])
            included_xsmiles = list(working_df["BB XSMILES"])
            included_xsmiles_valid = [
                bool_to_string_map[x] for x in list(working_df["XSMILES Valid"])
            ]
            for i in range(len(working_df)):
                data.append(
                    (
                        self.lib_name,
                        "A",  # TODO link this back to the the actual cycle that is encoded (it may not be "A")
                        included_xbbids[i],
                        included_xsmiles[i],
                        "",
                        included_xsmiles_valid[i],
                        "N",
                        "N",
                    )
                )

            df = pd.DataFrame(
                data,
                columns=[
                    "Library",
                    "Cycle",
                    "XBBID",
                    "XSMILES",
                    "Instruction",
                    "XSMILES Valid",
                    "Excluded",
                    "Exclude From Processing",
                ],
            )
            df.sort_values(
                by=[
                    "XSMILES Valid",
                    "Exclude From Processing",
                    "Excluded",
                ],
                ascending=[True, False, False],
                inplace=True,
            )

            df.to_csv(
                f"{xsmiles_directory}{self.lib_name}_A_BB_XSMILES.csv", # TODO link this back to the the actual cycle that is encoded (it may not be "A")
                index=False,
            )

            placeholder_df = pd.DataFrame(
                [
                    (
                        self.lib_name,
                        "A", # TODO link this back to the the actual cycle that is encoded (it may not be "A")
                        self.linker_cycle.default_smiles,
                        self.linker_cycle.placeholder_xsmiles_precursor,
                        self.linker_cycle.placeholder_xsmiles,
                    )
                ],
                columns=[
                    "Library",
                    "Cycle",
                    "Placeholder SMILES",
                    "Placeholder XSMILES Precursor",
                    "Placeholder XSMILES",
                ],
            )
            placeholder_df.to_csv(
                f"{xsmiles_directory}{self.lib_name}_A_placeholder_XSMILES.csv", # TODO link this back to the the actual cycle that is encoded (it may not be "A")
                index=False,
            )

        for cycle in self.cycles:
            data = []

            working_df = cycle.working_df.reset_index(drop=True)
            included_xbbids = list(working_df["XBBID"])
            included_xsmiles = list(working_df["BB XSMILES"])
            included_xsmiles_valid = [
                bool_to_string_map[x] for x in list(working_df["XSMILES Valid"])
            ]
            for i in range(len(working_df)):
                data.append(
                    (
                        self.lib_name,
                        cycle.cycle_name,
                        included_xbbids[i],
                        included_xsmiles[i],
                        "",
                        included_xsmiles_valid[i],
                        "N",
                        "N",
                    )
                )

            filtered_OUT_df_excluded_from_processing = deepcopy(cycle.filtered_OUT_df)
            filtered_OUT_df_excluded_from_processing = (
                filtered_OUT_df_excluded_from_processing[
                    filtered_OUT_df_excluded_from_processing[
                        "EXCLUDED_Excluded By User (XBBID)"
                    ]
                    == "Yes"
                ]
            )
            filtered_OUT_df_excluded_from_processing = (
                filtered_OUT_df_excluded_from_processing.reset_index(drop=True)
            )
            excluded_xbbids_excluded_from_processing = list(
                filtered_OUT_df_excluded_from_processing["EXCLUDED_XBBID"]
            )
            excluded_xsmiles_excluded_from_processing = list(
                filtered_OUT_df_excluded_from_processing["EXCLUDED_BB XSMILES"]
            )
            excluded_xsmiles_valid_excluded_from_processing = [
                bool_to_string_map[x]
                for x in list(
                    filtered_OUT_df_excluded_from_processing["EXCLUDED_XSMILES Valid"]
                )
            ]
            for i in range(len(filtered_OUT_df_excluded_from_processing)):
                data.append(
                    (
                        self.lib_name,
                        cycle.cycle_name,
                        excluded_xbbids_excluded_from_processing[i],
                        excluded_xsmiles_excluded_from_processing[i],
                        "FAIL<>",
                        excluded_xsmiles_valid_excluded_from_processing[i],
                        "Y",
                        "Y",
                    )
                )

            filtered_OUT_df_not_excluded_from_processing = deepcopy(
                cycle.filtered_OUT_df
            )
            filtered_OUT_df_not_excluded_from_processing = (
                filtered_OUT_df_not_excluded_from_processing[
                    filtered_OUT_df_not_excluded_from_processing[
                        "EXCLUDED_Excluded By User (XBBID)"
                    ]
                    != "Yes"
                ]
            )
            filtered_OUT_df_not_excluded_from_processing = (
                filtered_OUT_df_not_excluded_from_processing.reset_index(drop=True)
            )
            excluded_xbbids_not_excluded_from_processing = list(
                filtered_OUT_df_not_excluded_from_processing["EXCLUDED_XBBID"]
            )
            excluded_xsmiles_not_excluded_from_processing = list(
                filtered_OUT_df_not_excluded_from_processing["EXCLUDED_BB XSMILES"]
            )
            excluded_xsmiles_valid_not_excluded_from_processing = [
                bool_to_string_map[x]
                for x in list(
                    filtered_OUT_df_not_excluded_from_processing[
                        "EXCLUDED_XSMILES Valid"
                    ]
                )
            ]
            for i in range(len(filtered_OUT_df_not_excluded_from_processing)):
                data.append(
                    (
                        self.lib_name,
                        cycle.cycle_name,
                        excluded_xbbids_not_excluded_from_processing[i],
                        excluded_xsmiles_not_excluded_from_processing[i],
                        "",
                        excluded_xsmiles_valid_not_excluded_from_processing[i],
                        "Y",
                        "N",
                    )
                )
            df = pd.DataFrame(
                data,
                columns=[
                    "Library",
                    "Cycle",
                    "XBBID",
                    "XSMILES",
                    "Instruction",
                    "XSMILES Valid",
                    "Excluded",
                    "Exclude From Processing",
                ],
            )
            df.sort_values(
                by=[
                    "XSMILES Valid",
                    "Exclude From Processing",
                    "Excluded",
                ],
                ascending=[True, False, False],
                inplace=True,
            )

            df.to_csv(
                f"{xsmiles_directory}{self.lib_name}_{cycle.cycle_name}_BB_XSMILES.csv",
                index=False,
            )

            placeholder_df = pd.DataFrame(
                [
                    (
                        self.lib_name,
                        cycle.cycle_name,
                        cycle.reaction_branches[0].monosynthon_exemplar,
                        cycle.placeholder_xsmiles_precursor,
                        cycle.placeholder_xsmiles,
                    )
                ],
                columns=[
                    "Library",
                    "Cycle",
                    "Placeholder SMILES",
                    "Placeholder XSMILES Precursor",
                    "Placeholder XSMILES",
                ],
            )
            placeholder_df.to_csv(
                f"{xsmiles_directory}{self.lib_name}_{cycle.cycle_name}_placeholder_XSMILES.csv",
                index=False,
            )
        return None

    def enumerate_disynthons(
        self, disynthon_type: Tuple[str], enum_destination: str
    ) -> pa.Table:
        # Cycles should already be sorted in ascending order, but since this is critical, sort them again
        cycle_list = list(self.cycles)
        cycle_list.sort(key=lambda x: x.cycle_name)
        self.cycles = tuple(cycle_list)

        cycle_1 = disynthon_type[0]
        cycle_2 = disynthon_type[1]
        cycle_name_list = [cycle.cycle_name for cycle in self.cycles]
        # Initialize lists. these will be the same length as the number of cycles, and each element will be a 3-tuple
        with_instruction_rows = []
        no_instruction_rows = []

        expected_table_size = 1  # init
        for cycle in self.cycles:
            if cycle.cycle_name in disynthon_type:

                expected_table_size *= cycle.split_size
                with_instruction_xbbids = list(
                    cycle.xsmiles_dict_with_instruction.keys()
                )
                with_instruction_xsmiles_and_instruction = list(
                    cycle.xsmiles_dict_with_instruction.values()
                )
                no_instruction_xbbids = list(cycle.xsmiles_dict_no_instruction.keys())
                no_instruction_xsmiles_and_instruction = list(
                    cycle.xsmiles_dict_no_instruction.values()
                )
                # each row will be a 3-tuple of the form (XBBID, XSMILES, INSTRUCTION)
                with_instruction_rows.append(
                    [
                        tuple(
                            [
                                with_instruction_xbbids[i],
                                with_instruction_xsmiles_and_instruction[i][0],
                                with_instruction_xsmiles_and_instruction[i][1],
                            ]
                        )
                        for i in range(
                            len(list(cycle.xsmiles_dict_with_instruction.keys()))
                        )
                    ]
                )
                no_instruction_rows.append(
                    [
                        tuple(
                            [
                                no_instruction_xbbids[i],
                                no_instruction_xsmiles_and_instruction[i][0],
                                no_instruction_xsmiles_and_instruction[i][1],
                            ]
                        )
                        for i in range(
                            len(list(cycle.xsmiles_dict_no_instruction.keys()))
                        )
                    ]
                )
            else:
                with_instruction_rows.append([])
                no_instruction_rows.append([("Plc", cycle.placeholder_xsmiles, "")])

        with_instruction_combos = (
            generate_unique_combinations_with_at_least_one_instruction(
                with_instruction_rows=with_instruction_rows,
                no_instruction_rows=no_instruction_rows,
            )
        )
        no_instruction_combos = generate_combinations_with_no_instruction(
            no_instruction_rows=no_instruction_rows
        )

        if with_instruction_combos is not None:
            n_with_instruction = len(with_instruction_combos)
            combos = with_instruction_combos
            xsmiles_enum = parallel_map_concatenate_xsmiles(
                with_instruction_combos,
                cycle_name_list=cycle_name_list,
                use_instruction=True,
            )
            if no_instruction_combos is not None:
                n_without_instruction = len(no_instruction_combos)
                combos = combos + no_instruction_combos
                xsmiles_enum = xsmiles_enum + parallel_map_concatenate_xsmiles(
                    no_instruction_combos,
                    cycle_name_list=cycle_name_list,
                    use_instruction=False,
                )
            else:
                n_without_instruction = 0
        elif no_instruction_combos is not None:
            n_with_instruction = 0
            n_without_instruction = len(no_instruction_combos)
            combos = no_instruction_combos
            xsmiles_enum = pa.array(
                parallel_map_concatenate_xsmiles(
                    no_instruction_combos,
                    cycle_name_list=cycle_name_list,
                    use_instruction=False,
                )
            )
        else:
            raise Exception(
                "Both with_instruction_combos and no_instruction_combos are None"
            )

        # Make sure the number of combination matches what is expected based on the cycles' split sizes
        try:
            assert len(combos) == expected_table_size
        except:
            raise Exception(
                f"Unexpected number of combinations, {len(combos)}, were generated (expecting {expected_table_size})"
            )

        bb1_id_vector = [c[cycle_name_list.index(cycle_1)][0] for c in combos]
        bb2_id_vector = [c[cycle_name_list.index(cycle_2)][0] for c in combos]
        bb1_id_df = pd.DataFrame({"BBs": bb1_id_vector})
        bb2_id_df = pd.DataFrame({"BBs": bb2_id_vector})
        if enum_destination == "X-Viewer":
            instruction_vector = ["Y"] * n_with_instruction + [
                "N"
            ] * n_without_instruction
            bb1_xs_vector = pa.array(
                [c[cycle_name_list.index(cycle_1)][1] for c in combos]
            )
            bb2_xs_vector = pa.array(
                [c[cycle_name_list.index(cycle_2)][1] for c in combos]
            )
            bb1_instruction_vector = pa.array(
                [c[cycle_name_list.index(cycle_1)][2] for c in combos]
            )
            bb2_instruction_vector = pa.array(
                [c[cycle_name_list.index(cycle_2)][2] for c in combos]
            )
            all_bbs_df = query_candidate_bbs(xviewer_connection=self.xviewer_connection)
            bb1_smi_vector = pa.array(
                list(
                    bb1_id_df.join(all_bbs_df.set_index("XBBID"), on="BBs", how="left")[
                        "SMILES"
                    ]
                )
            )
            bb2_smi_vector = pa.array(
                list(
                    bb2_id_df.join(all_bbs_df.set_index("XBBID"), on="BBs", how="left")[
                        "SMILES"
                    ]
                )
            )
            names = [
                f"BB{disynthon_type[0]}",
                f"BB{disynthon_type[1]}",
                f"BB{disynthon_type[0]} SMILES",
                f"BB{disynthon_type[1]} SMILES",
                "SMILES",
                "Contains Instruction",
                f"BB{disynthon_type[0]} Instruction",
                f"BB{disynthon_type[1]} Instruction",
                f"BB{disynthon_type[0]} XS",
                f"BB{disynthon_type[1]} XS",
            ]
            tab = pa.Table.from_arrays(
                [
                    pa.array(bb1_id_vector),
                    pa.array(bb2_id_vector),
                    bb1_smi_vector,
                    bb2_smi_vector,
                    xsmiles_enum,
                    instruction_vector,
                    bb1_instruction_vector,
                    bb2_instruction_vector,
                    bb1_xs_vector,
                    bb2_xs_vector,
                ],
                names=names,
            )
        elif enum_destination == "S3":
            names = [
                "BB1",
                "BB2",
                "SMILES",
            ]
            tab = pa.Table.from_arrays(
                [
                    pa.array(bb1_id_vector),
                    pa.array(bb2_id_vector),
                    xsmiles_enum,
                ],
                names=names,
            )
        else:
            raise Exception(
                f"enum_destination must be either X-Viewer or S3, but {enum_destination} was input."
            )
        return tab

    def enumerate_and_write_instances(self, output_filepath: str, canonicalize_smiles: bool = False) -> None:
        # Cycles should already be sorted in ascending order, but since this is critical, sort them again
        cycle_list = list(self.cycles)
        cycle_list.sort(key=lambda x: x.cycle_name)
        self.cycles = tuple(cycle_list)

        cycle_name_list = [cycle.cycle_name for cycle in self.cycles]
        # Initialize lists. these will be the same length as the number of cycles, and each element will be a 3-tuple
        with_instruction_rows = []
        no_instruction_rows = []

        expected_table_size = 1  # init
        actual_table_size = 0  # init
        existing_data_behavior_for_output_parquet_files = "delete_matching"  # init
        k_for_parquet_file_basename_template = 1  # init

        for cycle in self.cycles:
            expected_table_size *= cycle.split_size
            with_instruction_xbbids = list(cycle.xsmiles_dict_with_instruction.keys())
            with_instruction_xsmiles_and_instruction = list(
                cycle.xsmiles_dict_with_instruction.values()
            )
            no_instruction_xbbids = list(cycle.xsmiles_dict_no_instruction.keys())
            no_instruction_xsmiles_and_instruction = list(
                cycle.xsmiles_dict_no_instruction.values()
            )
            # each row will be a 3-tuple of the form (XBBID, XSMILES, INSTRUCTION)
            with_instruction_rows.append(
                [
                    tuple(
                        [
                            with_instruction_xbbids[i],
                            with_instruction_xsmiles_and_instruction[i][0],
                            with_instruction_xsmiles_and_instruction[i][1],
                        ]
                    )
                    for i in range(
                        len(list(cycle.xsmiles_dict_with_instruction.keys()))
                    )
                ]
            )
            no_instruction_rows.append(
                [
                    tuple(
                        [
                            no_instruction_xbbids[i],
                            no_instruction_xsmiles_and_instruction[i][0],
                            no_instruction_xsmiles_and_instruction[i][1],
                        ]
                    )
                    for i in range(len(list(cycle.xsmiles_dict_no_instruction.keys())))
                ]
            )

        print(f"Expected size: {expected_table_size:,}")

        with_instruction_combos = (
            generate_unique_combinations_with_at_least_one_instruction(
                with_instruction_rows=with_instruction_rows,
                no_instruction_rows=no_instruction_rows,
            )
        )

        if with_instruction_combos is not None:
            xsmiles_enum = parallel_map_concatenate_xsmiles(
                    with_instruction_combos,
                    cycle_name_list=cycle_name_list,
                    use_instruction=True,
                )
            if canonicalize_smiles:
                xsmiles_enum = parallel_map_canonicalize_smiles_in_chunks(xsmiles_enum)

            xsmiles_enum = pa.array(xsmiles_enum)

            column_values = []
            column_names = []

            for j in range(len(self.cycles)):
                column_values.append(
                    pa.array([c[j][0] for c in with_instruction_combos])
                )
                column_names.append(f"BB{str(j+1)}")

            column_values = column_values + [xsmiles_enum]
            column_names = column_names + ["SMILES"]

            tab = pa.Table.from_arrays(
                column_values,
                names=column_names,
            )
            pda.write_dataset(
                tab,
                output_filepath,
                format="parquet",
                basename_template="part-"
                + str(k_for_parquet_file_basename_template)
                + "-{i}.parquet",
                max_rows_per_file=10000,
                max_rows_per_group=10000,
                max_open_files=1000,
                existing_data_behavior=existing_data_behavior_for_output_parquet_files,
            )
            print(
                f"Wrote {len(xsmiles_enum)} instances with an instruction to part-{k_for_parquet_file_basename_template}-X.parquet"
            )

            # Running total
            actual_table_size += len(xsmiles_enum)

            # Free up memory
            del tab
            del column_values
            del with_instruction_combos
            del xsmiles_enum

            # Next iteration, don't overwrite, just append
            existing_data_behavior_for_output_parquet_files = "overwrite_or_ignore"
            k_for_parquet_file_basename_template += 1

        else:
            pass

        largest_cycle_idx = np.argmax([len(ls) for ls in no_instruction_rows])
        largest_cycle_no_instruction_rows = no_instruction_rows[largest_cycle_idx]
        largest_cycle_no_instruction_size = len(largest_cycle_no_instruction_rows)

        # chunk the list of combinations
        n_chunks = min(
            max(math.ceil((expected_table_size - actual_table_size) / 10000000), 1),
            largest_cycle_no_instruction_size,
        )
        chunked_largest_cycle = np.array_split(
            np.array(largest_cycle_no_instruction_rows), n_chunks
        )  # List[List[Tuple]]

        for idx in range(n_chunks):
            no_instruction_rows_copy = deepcopy(no_instruction_rows)
            no_instruction_rows_copy[largest_cycle_idx] = chunked_largest_cycle[idx]
            no_instruction_combos = generate_combinations_with_no_instruction(
                no_instruction_rows=no_instruction_rows_copy
            )
            if no_instruction_combos is not None:
                xsmiles_enum = parallel_map_concatenate_xsmiles(
                        no_instruction_combos,
                        cycle_name_list=cycle_name_list,
                        use_instruction=False,
                    )
                if canonicalize_smiles:
                    xsmiles_enum = parallel_map_canonicalize_smiles_in_chunks(xsmiles_enum)
                xsmiles_enum = pa.array(xsmiles_enum)

                column_values = []
                column_names = []

                for j in range(len(self.cycles)):
                    column_values.append(
                        pa.array([c[j][0] for c in no_instruction_combos])
                    )
                    column_names.append(f"BB{str(j+1)}")

                column_values = column_values + [xsmiles_enum]
                column_names = column_names + ["SMILES"]

                tab = pa.Table.from_arrays(
                    column_values,
                    names=column_names,
                )
                pda.write_dataset(
                    tab,
                    output_filepath,
                    format="parquet",
                    basename_template="part-"
                    + str(k_for_parquet_file_basename_template)
                    + "-{i}.parquet",
                    max_rows_per_file=10000,
                    max_rows_per_group=10000,
                    max_open_files=1000,
                    existing_data_behavior=existing_data_behavior_for_output_parquet_files,
                )

                print(
                    f"Wrote chunk {idx+1} of {n_chunks} without an instruction ({len(xsmiles_enum)} instances) to part-{k_for_parquet_file_basename_template}-X.parquet"
                )

                # Running total
                actual_table_size += len(xsmiles_enum)

                # Free up memory
                del tab
                del column_values
                del no_instruction_combos
                del xsmiles_enum

                # Next iteration, don't overwrite, just append
                existing_data_behavior_for_output_parquet_files = "overwrite_or_ignore"
                k_for_parquet_file_basename_template += 1

            else:
                pass

        try:
            assert actual_table_size == expected_table_size
        except:
            raise Exception(
                f"Table size, {actual_table_size}, does not match expected table size, {expected_table_size}"
            )

        return None
