from typing import List, Dict, Union
from itertools import chain

from rdkit import Chem

from configs.inputs import FunctionalGroupConfig
from tools.parallel_map import parallel_map_generic


def has_substructure_match(mol: Chem.rdchem.Mol, smarts: str) -> bool:
    """
    Checks whether a molecule contains a substructure
    """

    if mol is None:
        return False
    else:
        return mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))


def parallel_map_has_substructure_match(
    mol_array: Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol], smarts: str
) -> List[bool]:
    """
    Parallel map implementation of has_substructure_match(mol, smarts).
    Consumes a list of molecules, and returns a list of bools, whether or
    not the substructure is present in each molecule.
    """

    if isinstance(mol_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(mol_array, has_substructure_match, smarts=smarts)
            )
        )
    else:
        return has_substructure_match(mol_array)


def count_substructure_matches(mol: Chem.rdchem.Mol, smarts: str) -> int:
    """
    Counts the number of times a substructure appears in a molecule
    """

    if mol is None:
        return 0
    else:
        return len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))


def parallel_map_count_substructure_matches(
    mol_array: Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol], smarts: str
) -> List[bool]:
    """
    Parallel map implementation of count_substructure_matches
    """

    if isinstance(mol_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(
                    mol_array, count_substructure_matches, smarts=smarts
                )
            )
        )
    else:
        return count_substructure_matches(mol_array, smarts)


def generate_fg_match_dict(
    fg: FunctionalGroupConfig, monosynthons: List[Chem.rdchem.Mol]
) -> Dict:
    """
    Generates a dictionary of the form {functional group name: mask}, where mask is a list of bools, corresponding to
    whether each monosynthon contains the functional group or not
    """
    return {fg.fg_name: parallel_map_has_substructure_match(monosynthons, fg.smarts)}


def parallel_map_generate_fg_match_dicts(
    fg_array: Union[List[FunctionalGroupConfig], FunctionalGroupConfig],
    monosynthons: List[Chem.rdchem.Mol],
) -> List[Dict]:
    """
    Parallel map implementation of generate_fg_match_dict
    """

    if isinstance(fg_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(
                    fg_array, generate_fg_match_dict, monosynthons=monosynthons
                )
            )
        )
    else:
        return generate_fg_match_dict(fg_array)
