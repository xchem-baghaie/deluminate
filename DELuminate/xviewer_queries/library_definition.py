from copy import copy

import pandas as pd

from configs.inputs import (
    ValidationStepConfig,
    ReactionBranchConfig,
    CycleConfig,
    LibraryConfig,
    LinkerCycleConfig,
)
from configs.library_definition import (
    SYMMETRICAL_MONO_BOC_TRIAMINE_LIBRARY_IDS,
    BOC_AMINE_SNAR_LIBRARY_IDS,
    LIBRARIES_TO_SKIP_PLACEHOLDER_VALIDATION,
    LINKER_DEFAULT_SMILES_FILEPATH,
)
from configs.xviewer import XViewer_Connection
from transforms.monosynthon_exemplar_builder import build_monosynthon_exemplar
from xviewer_queries.generic_query import generic_query
from xviewer_queries.reactions import (
    generate_dict_of_reactions,
    generate_dict_of_functional_groups,
)

def query_library_synonym(
    library_name: str,
    xviewer_connection:XViewer_Connection = None,
) -> str:
    if xviewer_connection is None:
        return None
    synonym_df = generic_query(
        query=f"""
        SELECT LIBRARY_NAME, LIBRARY_SYNONYM 
        FROM XCRX_LIBREG.VW_LIBRARIES
        WHERE LIBRARY_NAME = '{library_name}'
        """,
        connection=xviewer_connection,
    )
    if len(synonym_df) == 1:
        synonym = synonym_df["LIBRARY_SYNONYM"][0]
        if isinstance(synonym, str) and synonym != '':
            return synonym
    return None

def query_library_deck(
    library_name: str,
    xviewer_connection:XViewer_Connection = None,
) -> str:
    if xviewer_connection is None:
        return None
    deck_df = generic_query(
        query=f"""
        SELECT LIBRARY_NAME, LIBRARY_DECK 
        FROM XCRX_LIBREG.VW_LIBRARIES
        WHERE LIBRARY_NAME = '{library_name}'
        """,
        connection=xviewer_connection,
    )
    if len(deck_df) == 1:
        deck = deck_df["LIBRARY_DECK"][0]
        if isinstance(deck, str) and deck != '':
            return deck
    return None

def define_library(
    library_name: str,
    filename: str,
    compute_xsmiles: bool = False,
    is_presynthesis: bool = False,
    is_virtual: bool = False,
    use_xviewer_connection: bool = True,
) -> LibraryConfig:
    
    if use_xviewer_connection:
        xviewer_connection = XViewer_Connection()
        xviewer_connection.get_credentials()
    else:
        xviewer_connection = None

    #If we're using an XViewer connection, we can query the synonym. Otherwise the library will have no synonym.
    synonym = query_library_synonym(library_name, xviewer_connection=xviewer_connection)
    
    #If we're using an XViewer connection and this is postsynthesis, we can query the deck. Otherwise the deck will be looked up from the input file.
    deck = None
    if not is_presynthesis:
        deck = query_library_deck(library_name, xviewer_connection=xviewer_connection)
    
    if deck is None:
        lib_info_df = pd.read_excel(filename, sheet_name="Info").fillna("")
        deck = str(lib_info_df["Deck"][0])
        if deck == "":
            deck = None
    
    cycle_config_list = []
    functional_groups = generate_dict_of_functional_groups()

    # instantiate linker cycle
    if is_virtual:
        linker_df = pd.read_excel(f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{library_name}/postsynthesis/input/{library_name}_postsynthesis.xlsx", sheet_name="L").fillna("")
    else:
        linker_df = pd.read_excel(filename, sheet_name="L").fillna("")
    
    encoding_tag_set = linker_df["Encoding Tag Set"][0]

    linker_default_smiles_df = pd.read_csv(LINKER_DEFAULT_SMILES_FILEPATH).fillna("")
    linker_default_smiles_dict = {
        list(linker_default_smiles_df["XBBID"])[i]: list(
            linker_default_smiles_df["Default SMILES"]
        )[i]
        for i in range(len(linker_default_smiles_df))
    }

    xbbid_smiles_dict = {}
    for i in range(len(linker_df)):
        xbbid = linker_df["Linker XBBIDs"][i]
        if xbbid != "":
            override_smiles = linker_df["Override SMILES"][i]
            if override_smiles != "":
                xbbid_smiles_dict[xbbid] = override_smiles
            else:
                xbbid_smiles_dict[xbbid] = linker_default_smiles_dict[xbbid]

    # Set the default SMILES as the shortest string
    default_smiles = min(list(xbbid_smiles_dict.values()), key=len)

    assert default_smiles != ""

    linker_cycle = LinkerCycleConfig(
        encoding_tag_set,
        xbbid_smiles_dict,
        default_smiles,
        xviewer_connection=xviewer_connection,
    )

    for cycle in ["A", "B", "C", "D", "E"]:

        # Regenerated for each cycle because we manipulate the reactions in certain cases, so this ensures we start fresh each time
        reactions = generate_dict_of_reactions(
            compute_xsmiles=compute_xsmiles,
            is_presynthesis=is_presynthesis,
        )

        df = pd.read_excel(filename, sheet_name=cycle).fillna("")

        #For backwards compatibility
        df.rename(columns={"Reaction Branch":"Reaction Sequence ID"}, inplace=True)
        df.rename(columns={"Include BB In This Reaction Branch (Only One Branch per XBBID is Allowed)":"Include BB In This Reaction Sequence"},inplace=True)
        df.rename(columns={col_name:col_name.replace("Branches", "Sequences") for col_name in list(df.columns)}, inplace=True)
        df.rename(columns={col_name:col_name.replace("Branch", "Sequence") for col_name in list(df.columns)}, inplace=True)

        if library_name in SYMMETRICAL_MONO_BOC_TRIAMINE_LIBRARY_IDS:
            # Special case, requiring the removal of certain functional groups from the list of forbidden functional groups for particular reactions
            if cycle == "A":
                for reaction_name, reaction in reactions.items():
                    if reaction_name == "Reductive Amination":
                        setattr(
                            reaction,
                            "forbidden_fgs_reactant2",
                            [
                                fg
                                for fg in reaction.forbidden_fgs_reactant2
                                if fg.fg_name != "Poly-Amine"
                            ],
                        )
        
        if library_name in BOC_AMINE_SNAR_LIBRARY_IDS and cycle == "C":
            #Special case, where we want to ignore the fact that Boc diamines were used in an SNAr reaction even though some may deprotect in the reaction, leading to regioisomers
            for reaction_name, reaction in reactions.items():
                if reaction.rxn_name == "SNAr (Amines)":
                    setattr(
                        reaction,
                        "forbidden_fgs_reactant2",
                        [
                            fg
                            for fg in reaction.forbidden_fgs_reactant2
                            if fg.fg_name != "Boc Protected Amine"
                        ],
                    )            
        if len(df) > 0:
            # Diversity Selection options
            if "Perform Diversity Selection" in list(
                df.columns
            ) and "Target Number of BBs" in list(df.columns) and not is_virtual:
                perform_diversity_selection = df["Perform Diversity Selection"][0]
                if perform_diversity_selection == "Yes":
                    perform_diversity_selection = True
                elif perform_diversity_selection == "No":
                    perform_diversity_selection = False
                else:
                    raise Exception(
                        f"Unexpected value for perform_diversity selection, {perform_diversity_selection}"
                    )

                if perform_diversity_selection:
                    target_number_of_bbs = df["Target Number of BBs"][0]
                    try:
                        assert isinstance(target_number_of_bbs, int)
                    except:
                        raise Exception(
                            f"Target number of BBs, {target_number_of_bbs}, is not an integer"
                        )
                else:
                    target_number_of_bbs = None
            else:
                perform_diversity_selection = False
                target_number_of_bbs = None

            # Validation Sequence definition
            validation_df = df.loc[
                :,
                df.columns.intersection(
                    [
                        "Reaction Sequence Linked to Validation Sequence",
                        "Validation Sequence Step",
                        "Validation Assay Names (separated by commas, if multiple)",
                        "Validation Assay to Use When Results are Equal",
                        "Validation Cutoff",
                        "Include Unvalidated BBs for this Validation Sequence Step",
                    ]
                ),
            ]

            validation_df = validation_df.sort_values(
                [
                    "Reaction Sequence Linked to Validation Sequence",
                    "Validation Sequence Step",
                ],
                ascending=[True, True],
            ).reset_index(drop=True)
            validation_df = validation_df[
                validation_df["Reaction Sequence Linked to Validation Sequence"] != ""
            ]
            validation_df["Reaction Sequence Linked to Validation Sequence"] = [
                int(x)
                for x in list(
                    validation_df["Reaction Sequence Linked to Validation Sequence"]
                )
            ]
            validation_sequence_provided = len(validation_df) > 0

            # Whether to filter for the same deck as the library in this cycle
            string = (
                str(df["Filter by Library Deck in This Cycle"][0]).lower().capitalize()
            )

            if string == "Yes":
                filter_for_deck = True
            elif string == "No":
                filter_for_deck = False
            else:
                raise ValueError(
                    f'Unsupported value ({string}) was given in the "Filter by Library Deck in This Cycle" column for cycle {cycle}'
                )

            # Create dataframes for user defined FG inclusions or exclusions
            added_fg_df = df.loc[
                :,
                df.columns.intersection(
                    [
                        "Manually Added FGs (this is used for bi- or trifunctionals, enter each on a new line)",
                        "Add To Reaction Sequence (separated by commas, if multiple)",
                    ]
                ),
            ].reset_index(drop=True)
            excluded_fg_df = df.loc[
                :,
                df.columns.intersection(
                    [
                        "Manually Excluded FGs (enter each on a new line)",
                        "Exclude From Reaction Sequence (separated by commas, if multiple)",
                    ]
                ),
            ].reset_index(drop=True)

            # Reaction Sequence Definition
            reaction_df = df.loc[
                :,
                df.columns.intersection(
                    ["Reaction Sequence ID", "Reaction Sequence Step", "Reaction Name"]
                ),
            ]
            reaction_df = reaction_df.sort_values(
                ["Reaction Sequence ID", "Reaction Sequence Step"],
                ascending=[True, True],
            ).reset_index(drop=True)
            reaction_branch_list = []

            for branch in list(
                set(reaction_df["Reaction Sequence ID"])
            ):  # unique branch values
                if branch != "":
                    branch = int(branch)
                    branch_df = reaction_df[
                        reaction_df["Reaction Sequence ID"] == branch
                    ]
                    branch_df = branch_df.reset_index(drop=True)
                    manually_added_fgs = []
                    manually_excluded_fgs = []
                    for i in range(len(added_fg_df)):
                        if (
                            str(
                                added_fg_df[
                                    "Add To Reaction Sequence (separated by commas, if multiple)"
                                ][i]
                            ).count(str(branch))
                            == 1
                        ):
                            fg_name = added_fg_df[
                                "Manually Added FGs (this is used for bi- or trifunctionals, enter each on a new line)"
                            ][i]
                            if fg_name != "":
                                manually_added_fgs.append(functional_groups[fg_name])
                    for i in range(len(excluded_fg_df)):
                        if (
                            str(
                                excluded_fg_df[
                                    "Exclude From Reaction Sequence (separated by commas, if multiple)"
                                ][i]
                            ).count(str(branch))
                            == 1
                        ):
                            fg_name = excluded_fg_df[
                                "Manually Excluded FGs (enter each on a new line)"
                            ][i]
                            if fg_name != "":
                                manually_excluded_fgs.append(functional_groups[fg_name])

                    # get reaction_sequence
                    branch_step_list = []  # init
                    for step in list(
                        set(branch_df["Reaction Sequence Step"])
                    ):  # unique branch step values
                        if step != "":
                            step = int(step)
                            step_df = branch_df[
                                branch_df["Reaction Sequence Step"] == step
                            ].reset_index(drop=True)
                            if len(step_df) != 1:
                                raise ValueError(
                                    f"Check whether all reaction branches have UNIQUE reaction branch steps in cycle {cycle}. If they are not unique, it is recommended to split into multiple reaction branches"
                                )
                            reaction = copy(
                                reactions[step_df["Reaction Name"][0]]
                            )  # instantiate ReactionConfig object
                            reaction.__setattr__("branch_step", (branch, step))
                            if step == 1:
                                reaction.__setattr__(
                                    "allowed_fgs_reactant2",
                                    [
                                        fg
                                        for fg in reaction.allowed_fgs_reactant2
                                        if fg.fg_name
                                        not in [
                                            fg2.fg_name for fg2 in manually_excluded_fgs
                                        ]
                                    ],
                                )
                                reaction.__setattr__(
                                    "manually_added_fgs_reactant2", manually_added_fgs
                                )
                                reaction.__setattr__(
                                    "manually_excluded_fgs_reactant2",
                                    manually_excluded_fgs,
                                )
                            branch_step_list.append(reaction)
                    reaction_sequence = tuple(branch_step_list)

                    # get validation_sequence
                    if validation_sequence_provided:
                        validation_sequence_for_reaction_branch_df = validation_df[
                            validation_df[
                                "Reaction Sequence Linked to Validation Sequence"
                            ]
                            == branch
                        ].reset_index(drop=True)
                        validation_step_list = []
                        for validation_step in list(
                            set(
                                validation_sequence_for_reaction_branch_df[
                                    "Validation Sequence Step"
                                ]
                            )
                        ):
                            if validation_step != "":
                                validation_step = int(validation_step)
                                validation_step_df = (
                                    validation_sequence_for_reaction_branch_df[
                                        validation_sequence_for_reaction_branch_df[
                                            "Validation Sequence Step"
                                        ]
                                        == validation_step
                                    ].reset_index(drop=True)
                                )
                                if len(validation_step_df) != 1:
                                    raise ValueError(
                                        f"Check whether all validation branch steps have UNIQUE values in cycle {cycle}. If they are not unique, it is recommended to split into multiple validation steps"
                                    )
                                assay_name_str = validation_step_df[
                                    "Validation Assay Names (separated by commas, if multiple)"
                                ][0]
                                assay_names = [
                                    x.strip()
                                    for x in assay_name_str.split(",")
                                    if x.strip() != ""
                                ]
                                assay_to_use_when_equal = validation_step_df[
                                    "Validation Assay to Use When Results are Equal"
                                ][0]
                                validation_cutoff = float(
                                    validation_step_df["Validation Cutoff"][0]
                                )
                                string = (
                                    str(
                                        validation_step_df[
                                            "Include Unvalidated BBs for this Validation Sequence Step"
                                        ][0]
                                    )
                                    .lower()
                                    .capitalize()
                                )

                                if string == "Yes":
                                    include_unvalidated_bbs = True
                                elif string == "No":
                                    include_unvalidated_bbs = False
                                else:
                                    raise ValueError(
                                        f'Unsupported value ({string}) was given in the "Include Unvalidated BBs for this Validation Sequence Step" column for cycle {cycle}'
                                    )

                                validation_step_list.append(
                                    ValidationStepConfig(
                                        validation_step,
                                        assay_names,
                                        validation_cutoff,
                                        assay_to_use_when_equal,
                                        include_unvalidated_bbs,
                                        xviewer_connection=xviewer_connection,
                                    )
                                )
                        validation_sequence = tuple(validation_step_list)
                    else:
                        validation_sequence = None

                    if len(reaction_sequence[0].allowed_fgs_reactant2) > 0:
                        monosynthon_exemplar = build_monosynthon_exemplar(
                            reaction_sequence[0].allowed_fgs_reactant2[0],
                            reaction_sequence[0].manually_added_fgs_reactant2,
                        )  # use the first reaction, since that contains all of the info about the required FGs
                    else:
                        monosynthon_exemplar = ""

                    reaction_branch_list.append(
                        ReactionBranchConfig(
                            branch,
                            reaction_sequence,
                            validation_sequence,
                            monosynthon_exemplar,
                            xviewer_connection=xviewer_connection,
                        )
                    )

            reaction_branches = tuple(reaction_branch_list)

            # User-defined XBBID inclusions and exclusions:
            included_xbbid_df = (
                df.loc[
                    :,
                    df.columns.intersection(
                        [
                            "Include ONLY these XBBIDs",
                            "Include BB In This Reaction Sequence",
                        ]
                    ),
                ]
                .drop_duplicates(subset="Include ONLY these XBBIDs", keep="first")
                .reset_index(drop=True)
            )
            included_xbbid_df = included_xbbid_df[included_xbbid_df['Include ONLY these XBBIDs'] != '']
            included_xbbid_df = included_xbbid_df.reset_index(drop=True)
            excluded_xbbid_df = (
                df.loc[
                    :,
                    df.columns.intersection(
                        ["Exclude these XBBIDs from ALL Reaction Sequences"]
                    ),
                ]
                .drop_duplicates(
                    subset="Exclude these XBBIDs from ALL Reaction Sequences",
                    keep="first",
                )
                .reset_index(drop=True)
            )
            excluded_xbbid_df = excluded_xbbid_df[excluded_xbbid_df['Exclude these XBBIDs from ALL Reaction Sequences'] != '']
            excluded_xbbid_df = excluded_xbbid_df.reset_index(drop=True)

            forced_inclusion_bbs = list(included_xbbid_df["Include ONLY these XBBIDs"])
            forced_inclusion_bbs_reaction_branches = list(
                included_xbbid_df["Include BB In This Reaction Sequence"]
            )
            forced_exclusion_bbs = list(
                excluded_xbbid_df["Exclude these XBBIDs from ALL Reaction Sequences"]
            )

            if len(forced_inclusion_bbs) == 0 or (
                len(forced_inclusion_bbs) == 1 and forced_inclusion_bbs[0] == ""
            ):
                forced_inclusion_bbs = None
            if len(forced_exclusion_bbs) == 0 or (
                len(forced_exclusion_bbs) == 1 and forced_exclusion_bbs[0] == ""
            ):
                forced_exclusion_bbs = None

            if forced_inclusion_bbs is not None and forced_exclusion_bbs is not None:
                overlap_count = len(
                    [x for x in forced_inclusion_bbs if x in forced_exclusion_bbs]
                )
                if overlap_count > 0:
                    print(
                        f'Warning: {overlap_count} XBBIDs were present in both the "Include ONLY these XBBIDs" list and the "Exclude these XBBIDs from ALL Reaction Sequences" list!'
                    )

            if forced_inclusion_bbs is not None:
                forced_inclusion_bb_reaction_branch_map = {}
                for i in range(len(forced_inclusion_bbs)):
                    xbbid = forced_inclusion_bbs[i]
                    if "ound" in xbbid:
                        raise Exception(
                            f'"{xbbid}" was listed as a BB to include in cycle {cycle}!'
                        )
                    forced_inclusion_bb_reaction_branch_map[xbbid] = (
                        forced_inclusion_bbs_reaction_branches[i]
                    )
            else:
                forced_inclusion_bb_reaction_branch_map = None

            if forced_exclusion_bbs is not None:
                for i in range(len(forced_exclusion_bbs)):
                    xbbid = forced_exclusion_bbs[i]
                    if "ound" in xbbid:
                        raise Exception(
                            f'"{xbbid}" was listed as a BB to exclude from cycle {cycle}!'
                        )

            cycle_config_list.append(
                CycleConfig(
                    cycle,
                    library_name,
                    reaction_branches,
                    filter_for_deck=filter_for_deck,
                    forced_inclusion_bb_dict=forced_inclusion_bb_reaction_branch_map,
                    forced_exclusion_bbs=forced_exclusion_bbs,
                    perform_diversity_selection=perform_diversity_selection,
                    target_number_of_bbs=target_number_of_bbs,
                    xviewer_connection=xviewer_connection,
                )
            )

    cycles = tuple(cycle_config_list)

    library = LibraryConfig(
        lib_name = library_name,
        linker_cycle = linker_cycle,
        cycles = cycles,
        synonym = synonym,
        deck = deck,
        xviewer_connection=xviewer_connection,
    )

    if compute_xsmiles:
        library.compute_placeholder_xsmiles()
        if library_name not in LIBRARIES_TO_SKIP_PLACEHOLDER_VALIDATION:
            library.validate_placeholder_xsmiles()

    return library
