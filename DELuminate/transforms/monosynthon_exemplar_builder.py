from typing import List, Tuple, Union

from rdkit import Chem

from configs.inputs import FunctionalGroupConfig
from configs.monosynthon_exemplar import COMMON_OVERWRITES
from transforms.rdkit_mols import compute_mol_from_smiles, compute_smiles_from_mol
from transforms.react import react_bimolecular, react_unimolecular


def append_functional_group_leaving_no_handle(
    mols: Tuple[Chem.rdchem.Mol],
) -> Chem.rdchem.Mol:
    """
    Combines two RDKit molecules with CH2D groups directly, leaving no CH2D group in the final molecule
    """
    return react_bimolecular(mols, "[CX4:1][2H].[CX4:2][2H]>>[C:1][C:2]")


def append_functional_group_leaving_a_handle(
    mols: Tuple[Chem.rdchem.Mol],
) -> Chem.rdchem.Mol:
    """
    Combines two RDKit molecules with CH2D groups through a tertiary amine, leaving a CH2D group in the final molecule
    """
    return react_bimolecular(mols, "[CX4:1][2H].[CX4:2][2H]>>[C:1]CN(C[C:2])C[2H]")


def remove_benzylic_methyl_groups(mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    while mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4H3]c")):
        mol = react_unimolecular(mol, "[CX4H3][c:1]>>[c:1]")
    return mol


def build_monosynthon_exemplar(
    reactive_fg: FunctionalGroupConfig, additional_fgs: List[FunctionalGroupConfig]
) -> Union[str, None]:
    """
    Builds a monosynthon exemplar using a defined functional group,
    as well as any functional groups that were manually added by the user. It constructs these using the deuteromethyl (CH2D) group present on the example_smiles for each
    functional group, then removes any deuteromethyl groups from the final structure
    """

    example_smiles = reactive_fg.example_smiles
    mol = compute_mol_from_smiles(example_smiles)
    if len(additional_fgs) > 0:
        for i in range(len(additional_fgs)):
            mols = tuple(
                [mol, compute_mol_from_smiles(additional_fgs[i].example_smiles)]
            )
            if (
                i == len(additional_fgs) - 1
            ):  # If we're done after this reaction, leave no trace
                mol = append_functional_group_leaving_no_handle(mols)
            else:  # We still need to leave a CH2D group behind to react with the next manually added fg
                mol = append_functional_group_leaving_a_handle(mols)
    smiles = compute_smiles_from_mol(mol).replace("[2H]", "")
    mol = compute_mol_from_smiles(smiles)
    mol = remove_benzylic_methyl_groups(mol)
    out_smi = compute_smiles_from_mol(mol)
    if out_smi in list(COMMON_OVERWRITES.keys()):
        out_smi = compute_smiles_from_mol(
            compute_mol_from_smiles(COMMON_OVERWRITES[out_smi])
        )
    return out_smi
