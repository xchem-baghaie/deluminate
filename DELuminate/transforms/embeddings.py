from typing import List, Tuple, Union
from itertools import chain

import numpy as np
from tdc.chem_utils import MolConvert

from configs.visualizations import TDC_FEATURES
from tools.parallel_map import parallel_map_generic
from transforms.rdkit_mols import parallel_map_compute_mol_from_smiles

Num2DFeatures = dict(
    ECFP2=2048,
    ECFP4=2048,
    ECFP6=2048,
    MACCS=167,
    Daylight=2048,
    RDKit2D=200,  # Float valued inside range [0, 1]
    PubChem=881,
)


def transform_tdc_2d_features(
    smiles: Union[str, List[str]], type_2d_feature: str = "ECFP2"
) -> np.ndarray:
    """
    Args:
        smiles, str
        type_2d_feature, str

    Returns:
        np.ndarray: Numpy array of TDC fingerprints
    """

    if type_2d_feature is None or type_2d_feature == "":
        return smiles
    assert type_2d_feature in Num2DFeatures
    num_smiles = len(smiles) if isinstance(smiles, list) else 1
    converter = MolConvert(src="SMILES", dst=type_2d_feature)
    if type_2d_feature == "RDKit2D":
        return converter(smiles).reshape(num_smiles, -1)
    else:
        return np.array([int(x) for x in converter(smiles)]).reshape(num_smiles, -1)


def parallel_map_numpy_tdc_features(
    smiles_array: List[str], type_2d_feature: str = "ECFP2"
) -> np.ndarray:
    """
    Parallel map implementation of of transform_tdc_2d_features
    """

    feat_array = list(
        chain.from_iterable(
            parallel_map_generic(
                smiles_array,
                transform_tdc_2d_features,
                type_2d_feature=type_2d_feature,
            )
        )
    )
    return np.concatenate(feat_array, axis=0)


def parallel_map_get_embeddings(smiles_list: List[str], embedding: str) -> np.ndarray:
    """
    Generates fingerprints for a list of SMILES, using `embedding`.
    """

    if embedding in TDC_FEATURES:
        return parallel_map_numpy_tdc_features(smiles_list, type_2d_feature=embedding)
    else:
        mol_array = parallel_map_compute_mol_from_smiles(smiles_list)
        if embedding == "Morgan":
            from rdkit.Chem import AllChem

            morgan_bitvect_nested = parallel_map_generic(
                mol_array,
                AllChem.GetMorganFingerprintAsBitVect,
                radius=2,
            )
            return np.array(list(chain.from_iterable(morgan_bitvect_nested)))
        if embedding == "MHFP6":
            from mhfp.encoder import MHFPEncoder

            enc = MHFPEncoder()

            mol_nested = parallel_map_generic(mol_array, enc.encode_mol, min_radius=0)
            return np.array(list(chain.from_iterable(mol_nested)))  # , dtype=np.int32)
        elif embedding == "MAP4":
            from map4.map4 import MAP4Calculator

            enc = MAP4Calculator()

            mol_nested = parallel_map_generic(mol_array, np.vectorize(enc.calculate))
            return np.array(list(chain.from_iterable(mol_nested)))  # , dtype=np.int32)
        elif embedding == "MQN":
            from rdkit.Chem import rdMolDescriptors

            mol_nested = parallel_map_generic(mol_array, rdMolDescriptors.MQNs_)
            return np.array(list(chain.from_iterable(mol_nested)))
        elif embedding == "MXFP":
            return NotImplemented
        elif embedding == "ISIDA":
            return NotImplemented
        else:
            raise Exception(f"Unsupported embedding, {embedding}, was used")
