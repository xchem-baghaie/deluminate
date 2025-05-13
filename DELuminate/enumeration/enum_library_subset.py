from typing import List, Tuple, Dict
import random

import pandas as pd

from configs.inputs import CycleConfig, LinkerCycleConfig
from configs.library_definition import SYMMETRICAL_MONO_BOC_TRIAMINE_LIBRARY_IDS
from configs.xsmiles import ARTIFICIAL_ISOTOPES
from filters.substructures import filter_undesirable_substructures
from transforms.rdkit_mols import (
    parallel_map_compute_mol_from_smiles,
    parallel_map_compute_smiles_from_mol,
    compute_mol_from_smiles,
    compute_smiles_from_mol,
)
from transforms.react import (
    parallel_map_react_bimolecular,
    parallel_map_react_unimolecular,
    react_bimolecular,
)
from transforms.xsmiles import parallel_map_compute_xsmiles
from xviewer_queries.building_blocks import lookup_external_bbids
from xviewer_queries.functional_groups import generate_dict_of_functional_groups

def generate_combinations(
    df: pd.DataFrame,
    cycles: Tuple[CycleConfig],
    cycle_df_dict: Dict,
    num_instances: int,
) -> pd.DataFrame:
    """
    Generates `num_instances` combinations of building blocks such that the full list of building blocks for each cycle is exhausted before picking from it again.
    The building blocks are chosen in a different random order each time, to avoid situations where LCM(split_size_1, split_size_2) < `num_instances`, and the same
    building blocks from different cycles would always get paired up. The sampling method ensures that if `num_instances` >= max(split_sizes), every BB is used
    in at least one instance.
    """

    for cycle in cycles:
        k = 0
        cyc_df = cycle_df_dict[cycle.cycle_name]
        if len(cyc_df) == 0:
            raise Exception(f"Cycle {cycle.cycle_name} has no BBs! Exiting...")
        idxs = list(cyc_df.index.values)
        xbbids = []
        smiles = []
        xsmiles = []
        xsmiles_valid = []
        reaction_branch_ids = []
        reaction_steps = {}
        max_num_steps = max(
            [len(branch.reaction_sequence) for branch in cycle.reaction_branches]
        )
        for i in range(max_num_steps):
            reaction_steps[i + 1] = []
        while k < num_instances:
            random.seed(k)
            shuffled_idxs = random.sample(idxs, len(idxs))
            for i in shuffled_idxs:
                if k >= num_instances:
                    break
                else:
                    k += 1
                    xbbids.append(cyc_df["XBBID"][i])
                    smiles.append(cyc_df["SMILES"][i])
                    reaction_branch_ids.append(cyc_df["Reaction Sequence ID"][i])
                    if "BB XSMILES" in list(cyc_df.columns):
                        xsmiles.append(cyc_df["BB XSMILES"][i])
                        xsmiles_valid.append(cyc_df["XSMILES Valid"][i])
                    for j in range(max_num_steps):
                        reaction_steps[j + 1].append(
                            cyc_df["Reaction " + str(j + 1)][i]
                        )
        df["BB" + cycle.cycle_name + " ID"] = xbbids
        if cycle.xviewer_connection is not None:
            df["BB" + cycle.cycle_name + " EBBID"] = lookup_external_bbids(xbbids, cycle.xviewer_connection)
        df["BB" + cycle.cycle_name + " SMILES"] = smiles
        df["BB" + cycle.cycle_name + " Reaction Sequence"] = reaction_branch_ids
        if len(xsmiles) > 0:
            df["BB" + cycle.cycle_name + " XSMILES"] = xsmiles
            df["BB" + cycle.cycle_name + " XSMILES Valid"] = xsmiles_valid
        for i in range(max_num_steps):
            df["BB" + cycle.cycle_name + " Reaction " + str(i + 1)] = reaction_steps[
                i + 1
            ]
    return df


def filter_df_for_only_novel_instances(
    df: pd.DataFrame, cycle_names: List[str]
) -> pd.DataFrame:
    """
    Filters a dataframe for only instances which do not use XBBs in every cycle.
    """

    master_mask = []
    for i in range(len(df)):
        minor_mask = []
        for cycle_name in cycle_names:
            if "XBB" in df["BB" + cycle_name + " ID"][i]:
                minor_mask.append(True)
            else:
                minor_mask.append(False)
        master_mask.append(not all(minor_mask))
    df = df[master_mask]
    df = df.reset_index(drop=True)
    return df


def enumerate_library_instances(
    library_name: str,
    cycles: Tuple[CycleConfig],
    linker_cycle: LinkerCycleConfig,
    num_instances: int = 3000,
    qc_set: bool = False,
    is_virtual: bool = False,
    add_substructure_warnings: bool = False,
) -> pd.DataFrame:
    """
    Enumerates instances from the library.
    """

    # Sort cycles in ascending order
    cycle_list = list(cycles)
    cycle_list.sort(key=lambda x: x.cycle_name)
    cycles = tuple(cycle_list)

    cycle_df_dict = {}
    try:
        rand = 42 + int(
            library_name[len(library_name) - 3 :]
        )  # reproducible, but different for every library
    except:  # has a letter as the last character of the library name
        rand = 42 + int(library_name[:-1][len(library_name[:-1]) - 3 :])
    if qc_set:
        for cycle in cycles:
            rand += 29  # reproducible, but different for every cycle
            if cycle.forced_inclusion_bb_dict is not None:
                # We need to add in building blocks which were used in the library, but were filtered out by this script, so they are still tested in the QC comparison.
                additional_bbs = cycle.filtered_OUT_df[
                    cycle.filtered_OUT_df["XBBID"].isin(
                        list(cycle.forced_inclusion_bb_dict.keys())
                    )
                ]
                concat = pd.concat([cycle.working_df, additional_bbs])
                concat = concat.reset_index(drop=True)
                cycle_df_dict[cycle.cycle_name] = concat.sample(
                    frac=1, random_state=rand
                ).reset_index(drop=True)
            else:
                cycle_df_dict[cycle.cycle_name] = cycle.working_df.sample(
                    frac=1, random_state=rand
                ).reset_index(drop=True)
        num_instances = max(
            [len(bbs) for bbs in list(cycle_df_dict.values())]
        )  # reset to use all BBs at least once
    else:
        for cycle in cycles:
            rand += 29  # reproducible, but different for every cycle
            cycle_df_dict[cycle.cycle_name] = cycle.working_df.sample(
                frac=1, random_state=rand
            ).reset_index(drop=True)

    compute_xsmiles = "BB XSMILES" in list(cycle_df_dict[cycles[0].cycle_name].columns)

    df = pd.DataFrame(
        {
            "Library": [library_name] * num_instances,
        }
    ).fillna("")

    if compute_xsmiles:
        linker_df = linker_cycle.df[["XBBID", "LOCKSMILES", "BB XSMILES"]]
        if linker_cycle.encoding_tag_set is not None:
            if "A" not in [cycle.cycle_name for cycle in cycles]:
                linker_df.rename(
                    columns={
                        "XBBID": "BBA ID",
                        "LOCKSMILES": "BBA SMILES",
                        "BB XSMILES": "BBA XSMILES",
                    },
                    inplace=True,
                )
            else:
                raise Exception(
                    "The linker cycle is encoded, but cannot rename to cycle A because cycle A already exists in the library as a chemistry cycle."
                )
        else:
            linker_df.rename(
                columns={
                    "XBBID": "BBL ID",
                    "LOCKSMILES": "BBL SMILES",
                    "BB XSMILES": "BBL XSMILES",
                },
                inplace=True,
            )
    else:
        linker_df = linker_cycle.df[["XBBID", "LOCKSMILES"]]
        if linker_cycle.encoding_tag_set is not None:
            if "A" not in [cycle.cycle_name for cycle in cycles]:
                linker_df.rename(
                    columns={"XBBID": "BBA ID", "LOCKSMILES": "BBA SMILES"},
                    inplace=True,
                )
            else:
                raise Exception(
                    "The linker cycle is encoded, but cannot rename to cycle A because cycle A already exists in the library as a chemistry cycle."
                )
        else:
            linker_df.rename(
                columns={"XBBID": "BBL ID", "LOCKSMILES": "BBL SMILES"}, inplace=True
            )

    assert len(linker_df) > 0
    # Duplicate the linker df until we have enough rows...
    while len(linker_df) < num_instances:
        linker_df = pd.concat([linker_df, linker_df])
    # ...then downsample to `num_instances` rows
    linker_df = linker_df.sample(n=num_instances, random_state=rand).reset_index(
        drop=True
    )

    df = pd.concat([df, linker_df], axis=1)
    if linker_cycle.encoding_tag_set is not None:
        df["Mols"] = parallel_map_compute_mol_from_smiles(list(df["BBA SMILES"]))
    else:
        df["Mols"] = parallel_map_compute_mol_from_smiles(list(df["BBL SMILES"]))

    if is_virtual:
        working_df = generate_combinations(
            df=df,
            cycles=cycles,
            cycle_df_dict=cycle_df_dict,
            num_instances=num_instances,
        )

        working_df = filter_df_for_only_novel_instances(
            df=working_df,
            cycle_names=[cycle.cycle_name for cycle in cycles],
        )

        counter = 0

        done = len(working_df) >= num_instances
        while not done:
            counter += 1
            print(
                f"Generating new instances, round {counter}. Starting with {len(working_df)} instances"
            )
            updated_cycle_df_dict = {}
            for cycle in cycles:
                rand += 17  # reproducible, but different for every cycle
                updated_cycle_df_dict[cycle.cycle_name] = cycle.working_df.sample(
                    frac=1, random_state=rand
                ).reset_index(drop=True)

            new_instances_df = generate_combinations(
                df=df,
                cycles=cycles,
                cycle_df_dict=updated_cycle_df_dict,
                num_instances=num_instances,
            )

            new_instances_filtered_df = filter_df_for_only_novel_instances(
                df=new_instances_df,
                cycle_names=[cycle.cycle_name for cycle in cycles],
            )

            working_df = pd.concat([working_df, new_instances_filtered_df], axis=0)
            print(f"{len(working_df)} instances after concatenation")
            working_df.drop_duplicates(
                subset=[f"BB{cycle.cycle_name} ID" for cycle in cycles], inplace=True
            )
            print(f"{len(working_df)} instances after duplicate removal")

            if len(working_df) >= num_instances:
                done = True

            if counter > 5 and not done:
                out_path = "./instance_exemplars.csv"
                new_instances_df.to_csv(out_path, index=False)
                raise Exception(
                    f"Failed to generate novel instance exemplar combinations after 5 rounds. Writing to {out_path}"
                )

        if len(working_df) > num_instances:
            working_df = working_df.sample(
                num_instances, random_state=rand
            ).reset_index(drop=True)

        df = working_df

        assert len(df) == num_instances

    else:
        df = generate_combinations(
            df=df,
            cycles=cycles,
            cycle_df_dict=cycle_df_dict,
            num_instances=num_instances,
        )

    # Finish this for all cycles, then start looping again, so subsequent sorting does not interfere with the randomization
    intermediate_mols = {}
    for cycle in cycles:
        df = df.sort_values(
            by="BB" + cycle.cycle_name + " Reaction Sequence"
        ).reset_index(drop=True)
        df_list = []
        for reaction_branch in cycle.reaction_branches:
            branch_df = df[
                df["BB" + cycle.cycle_name + " Reaction Sequence"]
                == reaction_branch.reaction_branch_idx
            ].reset_index(drop=True)
            mols = list(branch_df["Mols"])
            bbs = parallel_map_compute_mol_from_smiles(
                list(branch_df["BB" + cycle.cycle_name + " SMILES"])
            )
            for i in range(len(reaction_branch.reaction_sequence)):
                rxn = reaction_branch.reaction_sequence[i]
                smarts = rxn.smarts.replace(
                    str(ARTIFICIAL_ISOTOPES["L"]), ""
                )  # If there are exotic isotopes in the smarts string, remove them.
                if compute_xsmiles:
                    for isotope_value in list(ARTIFICIAL_ISOTOPES.values()):
                        smarts = smarts.replace("[" + str(isotope_value), "[")
                if rxn.num_reactants == 2:
                    mols = parallel_map_react_bimolecular(mols, bbs, smarts)
                    for j in range(len(bbs)):
                        # If the reaction fails due to the incoming BB, fail the enumeration.
                        if (
                            compute_smiles_from_mol(
                                react_bimolecular(
                                    (
                                        compute_mol_from_smiles(rxn.example_reactant1),
                                        bbs[j],
                                    ),
                                    smarts,
                                    fail_if_no_reaction=True,
                                )
                            )
                            == ""
                        ):
                            mols[j] = compute_mol_from_smiles("")
                    if (
                        library_name in SYMMETRICAL_MONO_BOC_TRIAMINE_LIBRARY_IDS
                        and cycle.cycle_name == "A"
                    ):
                        # Workaround for diamine BB in cycle A with relative stereochemistry (XBB042192)
                        for j in range(len(bbs)):
                            if branch_df["BBA ID"][j] == "XBB042192":
                                mols[j] = compute_mol_from_smiles(
                                    "CC(C)(OC(N1C[C@H]2NCCN(C)[C@H]2C1)=O)C"
                                )
                else:
                    mols = parallel_map_react_unimolecular(mols, smarts)
            branch_df["Mols"] = mols
            df_list.append(branch_df)
        df = pd.concat(df_list).reset_index(drop=True)
        intermediate_mols[cycle.cycle_name] = mols
    instance_smiles_list = parallel_map_compute_smiles_from_mol(list(df["Mols"]))
    df["Instance SMILES"] = instance_smiles_list
    if add_substructure_warnings:
        fgs_to_search = [
            fg
            for fg in generate_dict_of_functional_groups().values()
            if fg.compound_proposal_rule in ["Warn", "Alert"]
        ]
        df = filter_undesirable_substructures(df, list(df["Mols"]), fg_list=fgs_to_search)
        df.drop(columns=["Include", "Exclude Due to Alert Substructure"],inplace=True)
    if compute_xsmiles:
        if linker_cycle.encoding_tag_set is not None:
            xsmiles_fragment_lists = [
                list(df["BBA XSMILES"])
            ]  # Will become a list of lists
        else:
            xsmiles_fragment_lists = [
                list(df["BBL XSMILES"])
            ]  # Will become a list of lists
        xsmiles_combo_list = []  # Will become a list of tuples
        xsmiles_valid_list = []  # Will become a list of bools
        for i in range(len(cycles)):
            if i == 0 and linker_cycle.encoding_tag_set is None:
                # If the linker cycle is not encoded, the first chemistry cycle XSMILES will have the linker portion embedded, but we already have the linker cycle XSMILES, so we ignore it in the first chemistry cycle
                xsmiles_fragment_lists.append(
                    [
                        xs.split(".")[1] + "."
                        for xs in list(df["BB" + cycles[i].cycle_name + " XSMILES"])
                    ]
                )
            else:
                xsmiles_fragment_lists.append(
                    list(df["BB" + cycles[i].cycle_name + " XSMILES"])
                )
        for i in range(len(df)):
            xsmiles_fragments = []
            for xsmiles_fragment_list in xsmiles_fragment_lists:
                xsmiles_fragments.append(xsmiles_fragment_list[i])
            xsmiles_combo_list.append(tuple(xsmiles_fragments))
            xsmiles_valid_for_each_cycle = []
            for cycle in cycles:
                xsmiles_valid_for_each_cycle.append(
                    list(df["BB" + cycle.cycle_name + " XSMILES Valid"])[i]
                )
            xsmiles_valid_list.append(all(xsmiles_valid_for_each_cycle))
        instance_xsmiles_list = parallel_map_compute_xsmiles(xsmiles_combo_list)
        df["Instance XSMILES"] = instance_xsmiles_list
        canonical_xsmiles_list = parallel_map_compute_smiles_from_mol(
            parallel_map_compute_mol_from_smiles(instance_xsmiles_list)
        )
        df["Instance Canonical XSMILES"] = canonical_xsmiles_list
        rxn_xsmiles_match = [
            instance_smiles_list[i] == canonical_xsmiles_list[i]
            and not instance_smiles_list[i] == ""
            for i in range(len(df))
        ]
        xsmiles_behavior_expected = [
            xsmiles_valid_list[i] == rxn_xsmiles_match[i] for i in range(len(df))
        ]
        df["RXN-XSMILES Match"] = rxn_xsmiles_match
        df["XSMILES Behavior Expected"] = xsmiles_behavior_expected

    instances = list(df["Mols"])
    df.drop(columns=["Mols"], inplace=True)
    total = len(df)
    try:
        assert total == num_instances
    except:
        raise Exception(
            f"Was expecting to enumerate {num_instances} instances, but only {total} were enumerated"
        )

    return df, instances
