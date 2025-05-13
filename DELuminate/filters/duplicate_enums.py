from typing import List

import pandas as pd
from rdkit import Chem

from transforms.rdkit_mols import (
    parallel_map_compute_smiles_from_mol,
    compute_mol_from_smiles,
)
from xviewer_queries.functional_groups import generate_dict_of_functional_groups

fg_dict = generate_dict_of_functional_groups()

carboxylic_acid_smarts_string = fg_dict["Carboxylic Acid"].smarts
acyclic_ester_smarts_string = fg_dict["Ester (Acyclic)"].smarts


def filter_duplicate_monosynthon_enumerations(
    df: pd.DataFrame, mols: List[Chem.rdchem.Mol]
) -> pd.DataFrame:
    """
    Labels rows of `df` which should not be kept because they enumerated the same monosynthon as another row.
    Special cases (halides, acids/esters) are accounted for, otherwise the first row will be kept and others will be labeled.
    Returns the labeled `df`.
    """

    df["Duplicate Monosynthon Enumeration Of"] = [""] * len(df)
    cansmis = parallel_map_compute_smiles_from_mol(mols)
    for i in range(len(df)):
        idx_i = df["XBBID"][i]
        bb_smiles_i = df["SMILES"][i]
        cansmi = cansmis[i]
        stereo = list(df["STEREOCHEMISTRYNAME"])[i]
        if cansmi != "":
            for j in range(len(df)):
                if i!=j:
                    idx_j = df["XBBID"][j]
                    if (
                        cansmis[j] == cansmi
                        and list(df["STEREOCHEMISTRYNAME"])[j] == stereo
                    ):
                        bb_smiles_j = df["SMILES"][j]
                        if (
                            bb_smiles_i.count("Cl") != cansmi.count("Cl")
                            or bb_smiles_i.count("Br") != cansmi.count("Br")
                            or bb_smiles_i.count("I") != cansmi.count("I")
                        ):  # halide is consumed during the reaction
                            if (
                                (
                                    bb_smiles_i.count("Cl") == 1
                                    and bb_smiles_j.count("Br") == 1
                                )
                                or (
                                    bb_smiles_i.count("Br") == 1
                                    and bb_smiles_j.count("I") == 1
                                )
                                or (
                                    bb_smiles_i.count("Cl") == 1
                                    and bb_smiles_j.count("I") == 1
                                )
                            ):  # Keep the more reactive halide (at index `j`)
                                if df["Include"][i] == "Yes" and df["Include"][j] == "Yes":
                                    df.loc[i, ["Include"]] = "No"
                                    df.loc[i, ["Duplicate Monosynthon Enumeration Of"]] = (
                                        idx_j
                                    )
                        if len(
                            compute_mol_from_smiles(bb_smiles_i).GetSubstructMatches(
                                Chem.MolFromSmarts(carboxylic_acid_smarts_string)
                            )
                        ) != len(
                            compute_mol_from_smiles(cansmi).GetSubstructMatches(
                                Chem.MolFromSmarts(carboxylic_acid_smarts_string)
                            )
                        ) or len(
                            compute_mol_from_smiles(bb_smiles_i).GetSubstructMatches(
                                Chem.MolFromSmarts(acyclic_ester_smarts_string)
                            )
                        ) != len(
                            compute_mol_from_smiles(cansmi).GetSubstructMatches(
                                Chem.MolFromSmarts(acyclic_ester_smarts_string)
                            )
                        ):  # carboxylate is consumed during the reaction
                            if (
                                len(
                                    compute_mol_from_smiles(bb_smiles_i).GetSubstructMatches(
                                        Chem.MolFromSmarts(acyclic_ester_smarts_string)
                                    )
                                )
                                == 1
                                and len(
                                    compute_mol_from_smiles(bb_smiles_j).GetSubstructMatches(
                                        Chem.MolFromSmarts(carboxylic_acid_smarts_string)
                                    )
                                )
                                == 0
                            ):  # Keep the carboxylic acid (at index `j`)
                                if df["Include"][i] == "Yes" and df["Include"][j] == "Yes":
                                    df.loc[i, ["Include"]] = "No"
                                    df.loc[i, ["Duplicate Monosynthon Enumeration Of"]] = (
                                        idx_j
                                    )
                        if (
                            j > i
                            and df["Include"][i] == "Yes"
                            and df["Include"][j] == "Yes"
                        ):  # Keep the BB which comes earlier in the list
                            df.loc[j, ["Include"]] = "No"
                            df.loc[j, ["Duplicate Monosynthon Enumeration Of"]] = idx_i
    return df
