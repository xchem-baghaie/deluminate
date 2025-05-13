from itertools import chain
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from tools.parallel_map import parallel_map_generic
from transforms.rdkit_mols import compute_smiles_from_mol

nonemol = Chem.MolFromSmiles("")

def react_unimolecular(
    mol: Chem.rdchem.Mol,
    smarts: str,
    fail_if_no_reaction: bool = False,
    react_until_all_sites_have_been_reacted: bool = False,
) -> Chem.rdchem.Mol:
    """
    Reacts one molecule

    Args:
        mol (Chem.rdchem.Mol), the reactant
        smarts (str), the SMARTS string of the reaction
        fail_if_no_reaction (bool), if True, will return an empty molecule if the reaction fails. Otherwise, it will return the starting material
        react_until_all_sites_have_been_reacted (bool), if True, will react until all sites have been reacted. Otherwise, it will only react once
    Returns:
        Chem.rdchem.Mol representation of the product. 
    """

    rxn = AllChem.ReactionFromSmarts(smarts)
    prods = rxn.RunReactants([mol])

    if len(prods) == 0:
        if fail_if_no_reaction:
            return nonemol
        else:
            return mol  # if the reaction fails, just return the starting material.
    
    if react_until_all_sites_have_been_reacted:
        while len(rxn.RunReactants([prods[0][0]])) > 0: #The first product can still react
            new_prods = [] #If there are already multiple products, we have to handle each one separately
            for i in range(len(prods)):
                new_prods = new_prods + [new_prod for new_prod in rxn.RunReactants([prods[0][0]])]
            prods = tuple(new_prods)
    
    if len(prods) == 1 or len(set([compute_smiles_from_mol(x[0]) for x in prods])) == 1:
        #Products are equivalent, so we can just return the first one
        prod = prods[0][0]
        Chem.SanitizeMol(prod)
        prod.UpdatePropertyCache()
        return prod
    else:
        return nonemol

def react_bimolecular(
    mols: Tuple[Chem.rdchem.Mol],
    smarts: str,
    fail_if_no_reaction: bool = False,
    react_until_all_sites_have_been_reacted: bool = False,
) -> Chem.rdchem.Mol:
    """
    Reacts two molecules

    Args:
        mols (Tuple[Chem.rdchem.Mol]): a two-tuple, first and second reactant
        smarts (str), the SMARTS string of the reaction
        fail_if_no_reaction (bool), if True, will return an empty molecule if the reaction fails. Otherwise, it will return the starting material
        react_until_all_sites_have_been_reacted (bool), if True, will react until all sites have been reacted. Otherwise, it will only react once
    Returns:
        Chem.rdchem.Mol representation of the product.
    """
    rxn = AllChem.ReactionFromSmarts(smarts)
    # If we are attempting to do a null reaction, the second reactant will be the nonemol. In this case, we just run the reaction on the reactant
    if compute_smiles_from_mol(mols[1]) != "":
        prods = rxn.RunReactants(mols)
    else:
        prods = rxn.RunReactants([mols[0]])
    if len(prods) == 0:
        if fail_if_no_reaction:
            return nonemol
        else:
            return mols[0]  # if reaction fails, just return the starting material.

    if react_until_all_sites_have_been_reacted:
        while len(rxn.RunReactants((prods[0][0], mols[1]))) > 0: #The first product can still react
            new_prods = [] #If there are already multiple products, we have to handle each one separately
            for i in range(len(prods)):
                new_prods = new_prods + [new_prod for new_prod in rxn.RunReactants((prods[i][0], mols[1]))]
            prods = tuple(new_prods)

    if len(prods) == 1 or len(set([compute_smiles_from_mol(x[0]) for x in prods])) == 1:
        #Products are equivalent, so we can just return the first one
        prod = prods[0][0]
        Chem.SanitizeMol(prod)
        prod.UpdatePropertyCache()
        return prod
    else:
        return nonemol  # if reaction is ambiguous, return the nonemol

def parallel_map_react_unimolecular(
    mol_array: List[Chem.rdchem.Mol],
    smarts: str,
    fail_if_no_reaction: bool = False,
    react_until_all_sites_have_been_reacted: bool = False,
) -> List[Chem.rdchem.Mol]:
    """
    Parallel map implementation of react_unimolecular(). Reacts a list of molecules in one reaction

    Args:
        mol_array (Chem.rdchem.Mol), the reactant
        smarts (str), the SMARTS string of the reaction
        fail_if_no_reaction (bool), if True, will return an empty molecule if the reaction fails. Otherwise, it will return the starting material
        react_until_all_sites_have_been_reacted (bool), if True, will react until all sites have been reacted. Otherwise, it will only react once
    Returns:
        A list of products, as Chem.rdchem.Mol
    """
    return list(
        chain.from_iterable(
            parallel_map_generic(mol_array,
                                 react_unimolecular,
                                 smarts=smarts,
                                 fail_if_no_reaction=fail_if_no_reaction,
                                 react_until_all_sites_have_been_reacted=react_until_all_sites_have_been_reacted
                                 )
        )
    )

def parallel_map_react_bimolecular(
    mol1_array: List[Chem.rdchem.Mol],
    mol2_array: List[Chem.rdchem.Mol],
    smarts: str,
    fail_if_no_reaction: bool = False,
    react_until_all_sites_have_been_reacted: bool = False,
) -> List[Chem.rdchem.Mol]:
    """
    Parallel map implementation of react_bimolecular(). Reacts two lists of molecules (row-wise) in one reaction

    Args:
        mol1_array (List[Chem.rdchem.Mol]), the list of first reactants
        mol1_array (List[Chem.rdchem.Mol]), the list of second reactants
        smarts (str), the SMARTS string of the reaction
        fail_if_no_reaction (bool), if True, will return an empty molecule if the reaction fails. Otherwise, it will return the starting material
        react_until_all_sites_have_been_reacted (bool), if True, will react until all sites have been reacted. Otherwise, it will only react once
    Returns:
        A list of products, as Chem.rdchem.Mol
    """
    try:
        assert len(mol1_array) == len(mol2_array)
    except AssertionError:
        raise AssertionError(
            f"Molecule arrays of differing lengths ({len(mol1_array)} and {len(mol2_array)}) were passed into parallel_map_react_bimolecular"
        )

    mol_array = [(mol1_array[i], mol2_array[i]) for i in range(len(mol1_array))]

    return list(
        chain.from_iterable(
            parallel_map_generic(mol_array,
                                 react_bimolecular,
                                 smarts=smarts,
                                 fail_if_no_reaction=fail_if_no_reaction,
                                 react_until_all_sites_have_been_reacted=react_until_all_sites_have_been_reacted
                                 )
        )
    )
