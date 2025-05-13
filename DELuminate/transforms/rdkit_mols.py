from functools import lru_cache
from itertools import chain
from typing import List, Union

from rdkit import Chem

from tools.parallel_map import parallel_map_generic

MAX_CACHE_SIZE = 256


@lru_cache(maxsize=MAX_CACHE_SIZE)
def compute_mol_from_smiles(smiles: str) -> Chem.rdchem.Mol:
    """
    Converts smile string to a Chem.rdchem.Mol description, keeping a cache of strings to reduce
    computation time
    Args:
        smiles (str): SMILES string, i.e. COC, which is to be converted
    Returns:
    Chem.rdchem.Mol representation of SMILES string
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # raise ValueError("Invalid SMILES string: %s" % smiles)
        mol = Chem.MolFromSmiles("")
    return mol


def parallel_map_compute_mol_from_smiles(
    smiles_array: Union[List[str], str],
) -> List[Chem.rdchem.Mol]:
    if isinstance(smiles_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(smiles_array, compute_mol_from_smiles)
            )
        )
    else:
        return compute_mol_from_smiles(smiles_array)


@lru_cache(maxsize=MAX_CACHE_SIZE)
def compute_smiles_from_mol(mol: Chem.rdchem.Mol) -> str:
    """
    Converts a Chem.rdchem.Mol to a SMILES string, keeping a cache of mols to reduce computation time

    Args:
        mol (Chem.rdchem.Mol): The molecule to be converted to a SMILES string
    Returns:
        str, the SMILES of the input molecule
    """
    if mol is None:
        return ""
    else:
        return Chem.MolToSmiles(mol).replace("\\\\", "\\")


def parallel_map_compute_smiles_from_mol(
    mol_array: Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol],
) -> List[str]:
    if isinstance(mol_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(mol_array, compute_smiles_from_mol)
            )
        )
    else:
        return compute_smiles_from_mol(mol_array)

def parallel_map_canonicalize_smiles_in_chunks(smiles_list: List[str], chunk_size:int = 100000) -> List[str]:
    canonical_smiles = []
    #Divide list into chunks of size N
    if len(smiles_list) > chunk_size:
        chunks = [smiles_list[i*chunk_size:(i+1)*chunk_size] for i in range((len(smiles_list)+chunk_size-1)//chunk_size)]
    else:
        chunks = [smiles_list]
    for chunk in chunks:
        canonical_smiles = canonical_smiles + parallel_map_compute_smiles_from_mol(
            parallel_map_compute_mol_from_smiles(chunk)
        )
    return canonical_smiles