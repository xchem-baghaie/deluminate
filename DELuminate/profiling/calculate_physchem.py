from typing import List, Dict
from itertools import chain

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors3D
from rdkit.Chem import AllChem

from configs.property_profiles import (
    RDKIT_PROPERTIES,
    PMI_XY_COLUMN_NAMES,
    PMI_RDS_COLUMN_NAMES,
    VALID_PROPERTIES,
)
from tools.parallel_map import parallel_map_generic
from transforms.rdkit_mols import (
    compute_smiles_from_mol,
    parallel_map_compute_mol_from_smiles,
)


def calculate_rdkit_properties(
    mol: Chem.rdchem.Mol, properties: List[str]
) -> np.ndarray:
    """
    Consumes a RDKit mol as input and returns a numpy array of computed properties from RDKit
    """

    descriptors = rdMolDescriptors.Properties(properties)
    if compute_smiles_from_mol(mol) != "":
        computed_properties = list(descriptors.ComputeProperties(mol))
        return computed_properties
    else:
        return [""] * (len(properties))


def parallel_map_calculate_rdkit_properties(
    mol_array: List[Chem.rdchem.Mol], properties: List[str] = RDKIT_PROPERTIES
) -> pd.DataFrame:
    """
    Parallel map implementation of calculate_rdkit_properties

    Consumes a list of n RDKit molecules, calculates m properties from RDKit, and returns an n x m pd.DataFrame of the properties
    """

    return pd.DataFrame(
        list(
            chain.from_iterable(
                parallel_map_generic(
                    mol_array, calculate_rdkit_properties, properties=properties
                )
            )
        ),
        columns=properties,
    )


def calculate_pmi(mol: Chem.rdchem.Mol) -> List:
    """
    Takes a RDKit mol as input, and returns a 10-element list, corresponding to 'PMI_X_Avg', 'PMI_Y_Avg', 'PMI_X_All', and 'PMI_Y_All', Rod_Avg, Disc_Avg, Sphere_Avg, Rod_All, Disc_All, Sphere_All for that molecule.
    """

    if compute_smiles_from_mol(mol) != "":
        mol = AllChem.AddHs(mol)
        confIds = AllChem.EmbedMultipleConfs(mol)
        multiconf_pmi1 = [
            round(Descriptors3D.NPR1(mol, confId), 5) for confId in confIds
        ]
        multiconf_pmi2 = [
            round(Descriptors3D.NPR2(mol, confId), 5) for confId in confIds
        ]
        assert len(multiconf_pmi1) == len(multiconf_pmi2)
        pmi1_avg = round(np.mean(multiconf_pmi1), 5)
        pmi2_avg = round(np.mean(multiconf_pmi2), 5)
        rod_avg = pmi2_avg - pmi1_avg
        disc_avg = 2 - (2 * pmi2_avg)
        sphere_avg = pmi1_avg + pmi2_avg - 1
        rod_all = [
            round(multiconf_pmi2[i] - multiconf_pmi1[i], 5)
            for i in range(len(multiconf_pmi1))
        ]
        disc_all = [
            round(2 - (2 * multiconf_pmi2[i]), 5) for i in range(len(multiconf_pmi1))
        ]
        sphere_all = [
            round(multiconf_pmi1[i] + multiconf_pmi2[i] - 1, 5)
            for i in range(len(multiconf_pmi1))
        ]
        return [
            pmi1_avg,
            pmi2_avg,
            multiconf_pmi1,
            multiconf_pmi2,
            rod_avg,
            disc_avg,
            sphere_avg,
            rod_all,
            disc_all,
            sphere_all,
        ]
    else:
        return [""] * 10


def parallel_map_calculate_pmi(mol_array: List[Chem.rdchem.Mol]) -> pd.DataFrame:
    """
    Parallel map implementation of calculate_pmi.

    Consumes a list of n RDKit molecules, and returns a n x 10 pd.DataFrame
    """

    return pd.DataFrame(
        list(chain.from_iterable(parallel_map_generic(mol_array, calculate_pmi))),
        columns=PMI_XY_COLUMN_NAMES + PMI_RDS_COLUMN_NAMES,
    )


def calculate_properties(
    df: pd.DataFrame,
    properties: Dict = VALID_PROPERTIES,
    mols: List[Chem.rdchem.Mol] = None,
) -> pd.DataFrame:
    """
    Calculates physchem properties on `mols`.
    Returns `df` with physchem property columns appended
    """

    if mols is None:
        mols = parallel_map_compute_mol_from_smiles(list(df["Instance SMILES"]))

    # Calculate any RDKit properties which the dataframe is missing
    properties_present = []
    for col in list(df.columns):
        if col in list(properties.keys()):
            val = properties[col]
            if not val in properties_present:
                properties_present.append(val)
    properties_to_calculate = [
        x for x in RDKIT_PROPERTIES if properties[x] not in properties_present
    ]
    if len(properties_to_calculate) > 0:
        property_df = parallel_map_calculate_rdkit_properties(
            mols, properties_to_calculate
        )
        df = pd.concat([df, property_df], axis=1)

    # Calculate PMI coordinates if the dataframe is missing any of them
    if not len(
        list(set(PMI_XY_COLUMN_NAMES + PMI_RDS_COLUMN_NAMES) & set(list(df.columns)))
    ) == len(PMI_XY_COLUMN_NAMES + PMI_RDS_COLUMN_NAMES):
        pmi_df = parallel_map_calculate_pmi(mols)
        df = pd.concat([df, pmi_df], axis=1)

    return df
