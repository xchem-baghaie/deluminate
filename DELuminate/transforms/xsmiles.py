from itertools import chain
from typing import Dict, List, Tuple, Union

from rdkit import Chem

from configs.xsmiles import (
    ARTIFICIAL_ISOTOPES,
    ARTIFICIAL_ELEMENTS,
    BLANK_PLACEHOLDER,
    CONNECTIVITY_HANDLES,
    MAX_NUMBER_OF_DISCONNECTIONS,
    TERMINATE_CHARACTER,
)
from tools.parallel_map import parallel_map_generic
from transforms.rdkit_mols import compute_smiles_from_mol, compute_mol_from_smiles


def replace_smiles_with_artificial_elements(
    smiles: str,
    artificial_elements: str,
    max_number_of_disconnections: int = MAX_NUMBER_OF_DISCONNECTIONS,
) -> str:
    for i in range(max_number_of_disconnections):
        smiles = smiles.replace(f"[{i}*]", "*")
    return smiles.replace("*", artificial_elements)


def set_isotopes(
    mol: Chem.rdchem.Mol,
    isotope_value: int,
    isotope_values_to_ignore: List[int] = list(ARTIFICIAL_ISOTOPES.values()),
) -> Chem.rdchem.Mol:
    """
    Sets isotopes of all atoms which do not have isotopes in isotope_values_to_ignore to a given isotope_value.

    There is a workaround implemented for handling nitro reduction reactions, in which we search for NH2 groups bonded to a carbon, where
    the N isotope is not set, but the carbon isotope is set to a value in isotope_values_to_ignore. in those cases, we set the nitrogen
    isotope to the same as the carbon isotope.
    """

    for match in mol.GetSubstructMatches(
        Chem.MolFromSmarts("[N!$([NX3][A]=,#[A$([!#6])])&$([NX3H2][#6])]")
    ):
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "N" and atom.GetIsotope() == 0:
                neighbors = atom.GetNeighbors()
                if len(neighbors) != 1:
                    raise Exception(
                        f"Primary amine had more than one neighbor in molecule with smiles {compute_smiles_from_mol(mol)}, which should not occur."
                    )
                atom.SetIsotope(neighbors[0].GetIsotope())

    for atom in mol.GetAtoms():
        if not atom.GetIsotope() in isotope_values_to_ignore:
            atom.SetIsotope(isotope_value)

    return mol


def reset_isotopes(mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """
    Sets isotopes of all atoms in a molecule to zero, which is the RDKit default.
    """

    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return mol


def parallel_map_set_isotopes(
    mol_array: Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol], isotope_value: int
) -> Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol]:
    if isinstance(mol_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(
                    mol_array,
                    set_isotopes,
                    isotope_value=isotope_value,
                )
            )
        )
    elif isinstance(mol_array, Chem.rdchem.Mol):
        return set_isotopes(mol_array)
    else:
        raise Exception(f"Unsupported type for mol_array ({type(mol_array)})")


def parallel_map_reset_isotopes(
    mol_array: Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol]
) -> Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol]:
    if isinstance(mol_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(
                    mol_array,
                    reset_isotopes,
                )
            )
        )
    elif isinstance(mol_array, Chem.rdchem.Mol):
        return reset_isotopes(mol_array)
    else:
        raise Exception(f"Unsupported type for mol_array ({type(mol_array)})")


def sort_bonds(bond) -> int:
    return (
        bond.GetBeginAtom().GetIsotope()
        + bond.GetBeginAtom().GetAtomicNum()
        + bond.GetEndAtom().GetIsotope()
        + bond.GetEndAtom().GetAtomicNum()
    )


def fragment_into_xsmiles_precursor(mol: Tuple, cycle_name: str) -> str:
    incompatible_elements = list(
        set([string[1:3] for string in list(CONNECTIVITY_HANDLES.keys())])
    )

    cycle_isotope = str(ARTIFICIAL_ISOTOPES[cycle_name])
    tracker = (
        []
    )  # This is used to track which cycle pairs have already been used, so we don't repeat the same dummy atom replacements for duplicate cycle pairs
    done = False
    isotope_list = [str(atom.GetIsotope()) for atom in mol.GetAtoms()]
    if cycle_isotope in isotope_list:
        counter = 0
        while not done:
            counter += 1
            if counter == 100:
                return TERMINATE_CHARACTER
            Chem.Kekulize(
                mol
            )  # Otherwise, if an aromatic bond is broken, the mol will be corrupted on the next iteration
            bond_list = list(mol.GetBonds())
            bond_list.sort(
                key=sort_bonds
            )  # so if there are multiple disconnections, they should be assigned in an order that is reproducible across cycles

            for bond in bond_list:
                isotope_1 = bond.GetBeginAtom().GetIsotope()
                isotope_2 = bond.GetEndAtom().GetIsotope()
                if isotope_1 != 0 and isotope_2 != 0 and isotope_1 != isotope_2:
                    cycle_1 = next(
                        k for k, v in ARTIFICIAL_ISOTOPES.items() if v == isotope_1
                    )  # inverse dict lookup
                    cycle_2 = next(
                        k for k, v in ARTIFICIAL_ISOTOPES.items() if v == isotope_2
                    )  # inverse dict lookup
                    cycles = [cycle_1, cycle_2]
                    cycles.sort()
                    cycle_pair = tuple(cycles)
                    mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
                    for a in mol.GetAtoms():
                        a.SetIsAromatic(False)
                    smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
                    fragment_list = smiles.split(".")
                    frags_to_keep = []
                    for fragment in fragment_list:
                        if fragment.count(cycle_isotope) > 0:
                            frags_to_keep.append(fragment)
                    if len(frags_to_keep) == 1:
                        smiles = frags_to_keep[0]
                    elif len(frags_to_keep) == 0:
                        print(
                            f"No fragments in {fragment_list} contained the expected isotope, {cycle_isotope}"
                        )
                        return TERMINATE_CHARACTER
                    else:
                        print(
                            f"More than one fragment in {fragment_list} contained the expected isotope, {cycle_isotope}"
                        )
                        return TERMINATE_CHARACTER
                    replaced_smiles = replace_smiles_with_artificial_elements(
                        smiles,
                        ARTIFICIAL_ELEMENTS[cycle_pair][tracker.count(cycle_pair)],
                    )
                    mol = compute_mol_from_smiles(replaced_smiles)
                    tracker.append(cycle_pair)
                    break
                elif (
                    bond_list.index(bond) == len(bond_list) - 1
                ):  # We exhausted all of the bonds in the molecule
                    done = True
        root_atom_idx = 0  # init
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in incompatible_elements:
                root_atom_idx = atom.GetIdx()
                break
        mol = reset_isotopes(mol)
        return replace_blank_placeholder(
            Chem.MolToSmiles(mol, rootedAtAtom=root_atom_idx, canonical=False)
        )
    else:
        return TERMINATE_CHARACTER


def parallel_map_fragment_into_xsmiles_precursors(
    mol_array: List[Chem.rdchem.Mol], cycle_name: str
) -> List[str]:
    return list(
        chain.from_iterable(
            parallel_map_generic(
                mol_array, fragment_into_xsmiles_precursor, cycle_name=cycle_name
            )
        )
    )


def replace_blank_placeholder(
    smiles: str,
    replacement_value: str = "[H]",
    blank_placeholder: str = BLANK_PLACEHOLDER,
) -> str:
    return smiles.replace(blank_placeholder, replacement_value)


def parallel_map_replace_blank_placeholder(
    smiles_array: List[str],
    replacement_value: str = "[H]",
    blank_placeholder: str = BLANK_PLACEHOLDER,
) -> List[str]:
    return list(
        chain.from_iterable(
            parallel_map_generic(
                smiles_array,
                replace_blank_placeholder,
                replacement_value=replacement_value,
                blank_placeholder=blank_placeholder,
            )
        )
    )


def add_connectivity_handles(
    smiles: str, connectivity_handles: Dict = CONNECTIVITY_HANDLES
):
    for element_pair, connectivity_handle in connectivity_handles.items():
        # Connectivity handles enclosed in parentheses do not yield RDKit-parsable structures
        # So, first, we replace those cases (with all possible bond types):
        smiles = smiles.replace(f"({element_pair})", connectivity_handle)
        smiles = smiles.replace(f"(-{element_pair})", f"-{connectivity_handle}")
        smiles = smiles.replace(f"(:{element_pair})", f":{connectivity_handle}")
        smiles = smiles.replace(f"(={element_pair})", f"={connectivity_handle}")
        smiles = smiles.replace(f"(#{element_pair})", f"#{connectivity_handle}")
        smiles = smiles.replace(f"(/{element_pair})", f"/{connectivity_handle}")
        smiles = smiles.replace(f"(\{element_pair})", f"\{connectivity_handle}")
        smiles = smiles.replace(f"(\\{element_pair})", f"\\{connectivity_handle}")
        smiles = smiles.replace(f"(\\\{element_pair})", f"\\\{connectivity_handle}")
        smiles = smiles.replace(f"(\\\\{element_pair})", f"\\\\{connectivity_handle}")
        # Finally, we replace all cases which are left:
        smiles = smiles.replace(element_pair, connectivity_handle)
    return smiles


def parallel_map_add_connectivity_handles(
    smiles_array: Union[List[str], str],
    connectivity_handles: Dict = CONNECTIVITY_HANDLES,
)-> Union[List[str], str]:
    if isinstance(smiles_array, List):
        return list(
            chain.from_iterable(
                parallel_map_generic(
                    smiles_array,
                    add_connectivity_handles,
                    connectivity_handles=connectivity_handles,
                )
            )
        )
    elif isinstance(smiles_array, Chem.rdchem.Mol):
        return add_connectivity_handles(smiles_array)
    else:
        raise Exception(f"Unsupported type for smiles_array ({type(smiles_array)})")


def compute_xsmiles(xsmiles_fragments: Tuple[str]) -> str:
    return "".join(xsmiles_fragments)


def parallel_map_compute_xsmiles(xsmiles_fragment_array: List[Tuple[str]]) -> List[str]:
    return list(
        chain.from_iterable(
            parallel_map_generic(
                xsmiles_fragment_array,
                compute_xsmiles,
            )
        )
    )


def parse_xsmiles_instructions(instructions: str) -> Tuple[List[str]]:
    """
    Parses an instruction string into individual instructions and their arguments.
    If an empty string is passed in, None is returned for both
    """
    if instructions != "":
        split_string = [string for string in instructions.split(">") if string != ""]
        instructions = [string.split("<")[0] for string in split_string]
        args = [string.split("<")[1].split("|") for string in split_string]
        try:
            assert len(instructions) == len(args)
        except:
            raise Exception(
                f"Different numbers of instructions ({len(instructions)}) and arguments ({len(args)}) were generated"
            )
        return instructions, args
    else:
        return None, None
