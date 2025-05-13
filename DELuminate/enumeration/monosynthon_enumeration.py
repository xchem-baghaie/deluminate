from typing import List, Tuple
from copy import deepcopy
import pandas as pd
from rdkit import Chem

from configs.inputs import ReactionConfig
from transforms.react import (
    parallel_map_react_unimolecular,
    parallel_map_react_bimolecular,
)
from transforms.rdkit_mols import (
    parallel_map_compute_smiles_from_mol,
)
from transforms.xsmiles import parallel_map_reset_isotopes, reset_isotopes
from tools.substructure_searching import parallel_map_has_substructure_match

nonemol = Chem.MolFromSmiles("")

def monosynthon_enum_cycle(
    cycle_name: str,
    reaction_sequence: Tuple[ReactionConfig],
    df: pd.DataFrame,
    bbs: List[Chem.rdchem.Mol],
    monosynthons: List[Chem.rdchem.Mol],
) -> Tuple[pd.DataFrame, List[Chem.rdchem.Mol]]:
    """
    Applies the same reaction sequence to all monosynthons using `bbs`. Returns the modified dataframe and monosynthons.
    """

    if "Failed Enumeration in Cycle" not in list(df.columns):
        df["Failed Enumeration in Cycle"] = [""] * len(df)
    if "Incompatible FG" not in list(df.columns):
        df["Incompatible FG"] = [""] * len(df)
    if "Capable of Reacting Multiple Times" not in list(df.columns):
        df["Capable of Reacting Multiple Times"] = [""] * len(df)

    for j in range(len(reaction_sequence)):
        rxn = reaction_sequence[j]
        for fg in rxn.forbidden_fgs_reactant2:
            mask = parallel_map_has_substructure_match(bbs, fg.smarts)
            for k in range(len(mask)):
                if mask[k]:
                    df.loc[k, ["Include"]] = "No"
                    df.loc[k, ["Incompatible FG"]] = (
                        f"{fg.fg_name} on building block (incompatible in cycle {cycle_name} {rxn.rxn_name})"
                    )
        for fg in rxn.forbidden_fgs_reactant1:
            mask = parallel_map_has_substructure_match(monosynthons, fg.smarts)
            for k in range(len(mask)):
                if mask[k]:
                    df.loc[k, ["Include"]] = "No"
                    df.loc[k, ["Incompatible FG"]] = (
                        f"{fg.fg_name} on DNA (incompatible in cycle {cycle_name} {rxn.rxn_name})"
                    )
        monosynthon_smiles_start = parallel_map_compute_smiles_from_mol(monosynthons)
        if rxn.num_reactants == 2:
            monosynthons = parallel_map_react_bimolecular(monosynthons, bbs, rxn.smarts)
            if rxn.intermediate_isolated:
                monosynthons_repeated = [""] * len(monosynthons)
            else:
                monosynthons_repeated = parallel_map_compute_smiles_from_mol(
                    parallel_map_react_bimolecular(monosynthons, bbs, rxn.smarts)
                )
        elif rxn.num_reactants == 1:
            monosynthons = parallel_map_react_unimolecular(monosynthons, rxn.smarts)
            if rxn.intermediate_isolated:
                monosynthons_repeated = [""] * len(monosynthons)
            else:
                monosynthons_repeated = parallel_map_compute_smiles_from_mol(
                    parallel_map_react_unimolecular(monosynthons, rxn.smarts)
                )
        else:
            raise ValueError(
                f"Reaction has an unsupported number of reactants, {rxn.num_reactants}"
            )
        monosynthon_smiles = parallel_map_compute_smiles_from_mol(monosynthons)

        for i in range(len(monosynthons)):
            if (
                "Blank" not in reaction_sequence[0].rxn_name
                and monosynthon_smiles[i] == monosynthon_smiles_start[i]
            ):
                df.loc[i, ["Include"]] = "No"
                if df["Failed Enumeration in Cycle"][i] == "":
                    df.loc[i, ["Failed Enumeration in Cycle"]] = cycle_name
            if not monosynthons_repeated[i] in ["", monosynthon_smiles[i]]:
                df.loc[i, ["Include"]] = "No"
                df.loc[i, ["Capable of Reacting Multiple Times"]] = "Yes"
    non_isotopic_monosynthons = deepcopy(monosynthons)
    df["SMILES after " + cycle_name] = parallel_map_compute_smiles_from_mol(
        parallel_map_reset_isotopes(non_isotopic_monosynthons)
    )
    for i in range(len(df)):
        if df["SMILES after " + cycle_name][i] == "":
            df.loc[i, ["Include"]] = "No"
            if df["Failed Enumeration in Cycle"][i] == "":
                df.loc[i, ["Failed Enumeration in Cycle"]] = cycle_name
    return df, monosynthons


def monosynthon_react_cycle(
    enum_cycle_unfiltered_df: pd.DataFrame,
    placeholder_cycle_name: str,
    placeholder_reaction_sequence: Tuple[ReactionConfig],
    placeholder_bb_exemplar: Chem.rdchem.Mol,
    monosynthons: List[Chem.rdchem.Mol],
) -> Tuple[pd.DataFrame, List[Chem.rdchem.Mol]]:
    """
    Applies the same reaction sequence to all monosynthons using the `placeholder_bb_exemplar`. Returns the modified dataframe and monosynthons.
    """

    df = enum_cycle_unfiltered_df

    if "Failed Enumeration in Cycle" not in list(df.columns):
        df["Failed Enumeration in Cycle"] = [""] * len(df)
    if "Incompatible FG" not in list(df.columns):
        df["Incompatible FG"] = [""] * len(df)
    if "Capable of Reacting Multiple Times" not in list(df.columns):
        df["Capable of Reacting Multiple Times"] = [""] * len(df)

    placeholder_bb_exemplar = reset_isotopes(placeholder_bb_exemplar)

    for j in range(len(placeholder_reaction_sequence)):
        rxn = placeholder_reaction_sequence[j]
        mask_list = []
        for fg in rxn.allowed_fgs_reactant1:
            mask_list.append(
                parallel_map_has_substructure_match(monosynthons, fg.smarts)
            )
        monosynthon_smiles_start = parallel_map_compute_smiles_from_mol(monosynthons)
        for k in range(len(monosynthons)):
            if monosynthon_smiles_start[k] != "" and not any(
                [mask_list[x][k] for x in range(len(mask_list))]
            ):
                df.loc[k, ["Include"]] = "No"
                if df["Incompatible FG"][k] == "":
                    df.loc[k, ["Incompatible FG"]] = (
                        f"No compatible functional group for cycle {placeholder_cycle_name} {rxn.rxn_name} was found (expecting one of {[fg.fg_name for fg in rxn.allowed_fgs_reactant1]})"
                    )

        for fg in rxn.forbidden_fgs_reactant1:
            mask = parallel_map_has_substructure_match(monosynthons, fg.smarts)
            for k in range(len(mask)):
                if mask[k]:
                    df.loc[k, ["Include"]] = "No"
                    df.loc[k, ["Incompatible FG"]] = (
                        f"{fg.fg_name} (incompatible in cycle {placeholder_cycle_name} {rxn.rxn_name})"
                    )
        if rxn.num_reactants == 2:
            monosynthons = parallel_map_react_bimolecular(
                monosynthons, [placeholder_bb_exemplar] * len(monosynthons), rxn.smarts
            )
            if rxn.intermediate_isolated:
                monosynthons_repeated = [""] * len(monosynthons)
            else:
                monosynthons_repeated = parallel_map_compute_smiles_from_mol(
                    parallel_map_react_bimolecular(
                        monosynthons,
                        [placeholder_bb_exemplar] * len(monosynthons),
                        rxn.smarts,
                    )
                )
        elif rxn.num_reactants == 1:
            monosynthons = parallel_map_react_unimolecular(monosynthons, rxn.smarts)
            if rxn.intermediate_isolated:
                monosynthons_repeated = [""] * len(monosynthons)
            else:
                monosynthons_repeated = parallel_map_compute_smiles_from_mol(
                    parallel_map_react_unimolecular(monosynthons, rxn.smarts)
                )
        else:
            raise ValueError(
                f"Reaction has an unsupported number of reactants, {rxn.num_reactants}"
            )
        monosynthon_smiles = parallel_map_compute_smiles_from_mol(monosynthons)
        for i in range(len(monosynthons_repeated)):
            if (
                "Blank" not in placeholder_reaction_sequence[0].rxn_name
                and monosynthon_smiles[i] == monosynthon_smiles_start[i]
            ):
                df.loc[i, ["Include"]] = "No"
                if df["Failed Enumeration in Cycle"][i] == "":
                    df.loc[i, ["Failed Enumeration in Cycle"]] = placeholder_cycle_name
            if not monosynthons_repeated[i] in ["", monosynthon_smiles[i]]:
                df.loc[i, ["Include"]] = "No"
                df.loc[i, ["Capable of Reacting Multiple Times"]] = "Yes"
    non_isotopic_monosynthons = deepcopy(monosynthons)
    df["SMILES after " + placeholder_cycle_name] = parallel_map_compute_smiles_from_mol(
        parallel_map_reset_isotopes(non_isotopic_monosynthons)
    )
    for i in range(len(df)):
        if df["SMILES after " + placeholder_cycle_name][i] == "":
            df.loc[i, ["Include"]] = "No"
            if df["Failed Enumeration in Cycle"][i] == "":
                df.loc[i, ["Failed Enumeration in Cycle"]] = placeholder_cycle_name
    return df, monosynthons
