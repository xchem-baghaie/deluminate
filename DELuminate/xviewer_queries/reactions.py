from typing import Dict

import pandas as pd

from configs.inputs import ReactionConfig
from configs.functional_groups import FGS_TO_IGNORE_FOR_PRESYNTHESIS, FGS_TO_IGNORE_FOR_POSTSYNTHESIS
from xviewer_queries.functional_groups import generate_dict_of_functional_groups


def generate_dict_of_reactions(
    filepath: str = "./DELuminate/xviewer_queries/Reactions.csv",
    compute_xsmiles: bool = False,
    is_presynthesis: bool = False,
) -> Dict:
    """
    Constructs a dictionary of the form {rxn_name: rxn object}, to mimic a table of reactions in X-Viewer.
    Currently, this pulls from a static file, which is generated manually.
    """
    df = pd.read_csv(filepath).fillna("")
    allowed_1 = []
    allowed_2 = []
    forbidden_1 = []
    forbidden_2 = []
    product_fgs = []
    associated_validation_assays_reactant_1 = []
    associated_validation_assays_reactant_2 = []
    for i in range(len(df)):
        allowed_1.append(
            [x.strip() for x in df["ALLOWED_FGS_REACTANT_1"][i].split(",") if x != ""]
        )
        allowed_2.append(
            [x.strip() for x in df["ALLOWED_FGS_REACTANT_2"][i].split(",") if x != ""]
        )
        forbidden_1.append(
            [x.strip() for x in df["FORBIDDEN_FGS_REACTANT_1"][i].split(",") if x != ""]
        )
        forbidden_2.append(
            [x.strip() for x in df["FORBIDDEN_FGS_REACTANT_2"][i].split(",") if x != ""]
        )
        product_fgs.append(
            [x.strip() for x in df["PRODUCT_FGS"][i].split(",") if x != ""]
        )
        associated_validation_assays_reactant_1.append(
            [
                x.strip()
                for x in df["ASSOCIATED_VALIDATION_ASSAYS_REACTANT_1"][i].split(",")
                if x != ""
            ]
        )
        associated_validation_assays_reactant_2.append(
            [
                x.strip()
                for x in df["ASSOCIATED_VALIDATION_ASSAYS_REACTANT_2"][i].split(",")
                if x != ""
            ]
        )
    df["Allowed_FGs_Reactant_1"] = allowed_1
    df["Allowed_FGs_Reactant_2"] = allowed_2
    df["Forbidden_FGs_Reactant_1"] = forbidden_1
    df["Forbidden_FGs_Reactant_2"] = forbidden_2
    df["Product_FGs"] = product_fgs
    rxn_dict = {}
    fg_dict = generate_dict_of_functional_groups()
    for i in range(len(df)):
        rxn_name = df["REACTION_NAME"][i]
        smarts = df["REACTION_SMARTS"][i]
        if compute_xsmiles:
            from configs.xsmiles import ARTIFICIAL_ISOTOPES

            if smarts == "[NH1:5]([CX4H3:8])[CX3:6](=S)[N:7][#6:9].[NX3H1$([NX3H1][NX3H2]):1]([NX3H2:2])[CX3:3](=O)[#6:4]>>[C:8][N:5]1[C:3]([#6:4])=[N:2][N:1]=[C:6]([N:7][#6:9])1":
                # Workaround for XLIB0108 because triazole formation uses a reactant that is symmetric, except for its isotopes (used to generate XSMILES)
                # So, here, we specify which carbon is from the linker side of the reactant, to help distinguish between them.
                smarts = f'[NH1:5]([{ARTIFICIAL_ISOTOPES["L"]}CX4H3:8])[CX3:6](=S)[N:7][#6:9].[NX3H1$([NX3H1][NX3H2]):1]([NX3H2:2])[CX3:3](=O)[#6:4]>>[C:8][N:5]1[C:3]([#6:4])=[N:2][N:1]=[C:6]([N:7][#6:9])1'
            elif smarts == "[NH1:5]([CX4H3:8])[CX3:6](=S)[N:7][#6:9].[NX3H1$([NX3H1][NX3H2]):1]([NX3H2:2])[CX3:3](=O)[#6:4]>>[#6:9][N:7]1[C:3]([#6:4])=[N:2][N:1]=[C:6]([N:5][C:8])1":
                #Same as above, but for regioisomer 2
                smarts = f'[NH1:5]([{ARTIFICIAL_ISOTOPES["L"]}CX4H3:8])[CX3:6](=S)[N:7][#6:9].[NX3H1$([NX3H1][NX3H2]):1]([NX3H2:2])[CX3:3](=O)[#6:4]>>[#6:9][N:7]1[C:3]([#6:4])=[N:2][N:1]=[C:6]([N:5][C:8])1'
        allowed_fgs_reactant1 = [fg_dict[fg] for fg in df["Allowed_FGs_Reactant_1"][i]]
        if is_presynthesis:
            # If we are in presynthesis mode, certain functional groups should not be flagged, so remove them from the list here
            forbidden_fgs_reactant1 = [
                fg_dict[fg]
                for fg in df["Forbidden_FGs_Reactant_1"][i]
                if fg not in FGS_TO_IGNORE_FOR_PRESYNTHESIS
            ]
        else:
            # If we are not in presynthesis mode, other functional groups should not be flagged, so remove them from the list here
            forbidden_fgs_reactant1 = [
                fg_dict[fg]
                for fg in df["Forbidden_FGs_Reactant_1"][i]
                if fg not in FGS_TO_IGNORE_FOR_POSTSYNTHESIS
            ]
        example_reactant1 = df["EXAMPLE_REACTANT_1"][i]
        if df["INTERMEDIATE_ISOLATED"][i] == "Y":
            intermediate_isolated = True
        elif df["INTERMEDIATE_ISOLATED"][i] == "N":
            intermediate_isolated = False
        else:
            raise Exception(
                "Unknown value in the INTERMEDIATE_ISOLATED column of the Reaction table"
            )
        allowed_fgs_reactant2 = [fg_dict[fg] for fg in df["Allowed_FGs_Reactant_2"][i]]
        if is_presynthesis:
            # If we are in presynthesis mode, certain functional groups should not be flagged, so remove them from the list here
            forbidden_fgs_reactant2 = [
                fg_dict[fg]
                for fg in df["Forbidden_FGs_Reactant_2"][i]
                if fg not in FGS_TO_IGNORE_FOR_PRESYNTHESIS
            ]
        else:
            # If we are not in presynthesis mode, other functional groups should not be flagged, so remove them from the list here
            forbidden_fgs_reactant2 = [
                fg_dict[fg]
                for fg in df["Forbidden_FGs_Reactant_2"][i]
                if fg not in FGS_TO_IGNORE_FOR_POSTSYNTHESIS
            ]
        example_reactant2 = df["EXAMPLE_REACTANT_2"][i]
        if example_reactant2 == "":
            example_reactant2 = None
        example_product = df["EXAMPLE_PRODUCT"][i]
        if example_product == "":
            example_product = None
        product_fgs = [fg_dict[fg] for fg in df["Product_FGs"][i]]
        rxn = ReactionConfig(
            rxn_name = rxn_name,
            smarts = smarts,
            allowed_fgs_reactant1 = allowed_fgs_reactant1,
            forbidden_fgs_reactant1 = forbidden_fgs_reactant1,
            example_reactant1 = example_reactant1,
            intermediate_isolated = intermediate_isolated,
            associated_validation_assays_reactant_1= associated_validation_assays_reactant_1,
            associated_validation_assays_reactant_2 = associated_validation_assays_reactant_2,
            allowed_fgs_reactant2 = allowed_fgs_reactant2,
            forbidden_fgs_reactant2 =forbidden_fgs_reactant2,
            product_fgs = product_fgs,
            example_reactant2 = example_reactant2,
            example_product = example_product,
        )
        rxn_dict[rxn_name] = rxn
    return rxn_dict

if __name__ == '__main__':
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np
    from time import time

    N = 100
    BATCH_SIZE = 10000
    rxn_dict = generate_dict_of_reactions()
    rxn_dict = {k:v for k, v in rxn_dict.items() if v.example_reactant2 is not None and k not in [
        'SNAr onto pyrimidine (4-Cl displacement) (Pyrimidine BBs)',
        'XLIB0095 SNAr of Dichloropyrimidine Core',
        'Michael Addition (Vinyl Sulfone BBs)',
        'Michael Addition (Vinyl Sulfone BBs) (Aliphatic Amine Only)',
        ] and k[0].lower() != 'z'}
    out = {
        'Reaction Name':[],
        'Enumeration Method':[],
        'Mean Processing Time (ms)':[],
        'Std Processing Time (ms)':[],
        'SMARTS':[],
        'SMARTS String Length':[],
        'Reactant 1 SMILES String Length':[],
        'Mean Reactant 2 SMILES String Length':[],
        'Std Reactant 2 SMILES String Length':[],
        'Reactant 1 HAC':[],
        'Mean Reactant 2 HAC':[],
        'Std Reactant 2 HAC':[],
        'N Batches':[],
        'Batch Size':[],
    }
    for rxn in list(rxn_dict.values()):
        smiles_1 = [rxn.example_reactant1]*BATCH_SIZE
        try:
            smiles_2_df = pd.read_csv(f'/mnt/p/discovery chemistry/people/ryan/tech dev - cheminformatics/enumeration/compatible_bbs_by_reaction/{rxn.rxn_name}.csv')
        except:
            try:
                smiles_2_df = pd.read_csv(f'/mnt/p/discovery chemistry/people/ryan/tech dev - cheminformatics/enumeration/compatible_bbs_by_reaction/{rxn.rxn_name} (Primary Amines).csv')
            except:
                print(f'{rxn.rxn_name} has no example SMILES data')
                continue
        smiles_2_ls = list(smiles_2_df['SMILES'])
        while len(smiles_2_ls) < BATCH_SIZE:
            smiles_2_ls = smiles_2_ls + smiles_2_ls
        smiles_2 = smiles_2_ls[:BATCH_SIZE]
        try:
            assert len(smiles_2) == len(smiles_1)
        except:
            raise Exception(f'{len(smiles_2)}, {len(smiles_1)}')
        smiles_2_lengths = [len(x) for x in smiles_2]

        rxn_times = []
        for _ in range(N):
            #Process Reaction
            t0 = time()
            mols_1 = [Chem.MolFromSmiles(s) for s in smiles_1]
            mols_2 = [Chem.MolFromSmiles(s) for s in smiles_2]
            reactor = AllChem.ReactionFromSmarts(rxn.smarts)
            rxn_prod_smis = []
            for i in range(BATCH_SIZE):
                prods = reactor.RunReactants((mols_1[i], mols_2[i]))
                try:
                    rxn_prod_smis.append(Chem.MolToSmiles(prods[0][0]))
                except:
                    rxn_prod_smis.append('')
            t1 = time()
            rxn_times.append((t1-t0)*1000)
            assert rxn_prod_smis.count('') <= len(rxn_prod_smis)/10 #No more than 10% failed to react
        out['Enumeration Method'].append('Reaction')
        out['Mean Processing Time (ms)'].append(round(np.mean(rxn_times),8))
        out['Std Processing Time (ms)'].append(round(np.std(rxn_times),8))

        xs_times = []
        for _ in range(N):
            #Process String Concatenation
            t2 = time()
            xs_prod_smis = []
            for i in range(BATCH_SIZE):
                try:
                    xs_prod_smis.append(''.join([smiles_1[i], smiles_2[i]]))
                except:
                    #This should never happen, but implemented this as try/except to be a fair comparison vs. reaction-based enumeration, where we do have to handle those situations
                    xs_prod_smis.append('')
            t3 = time()
            xs_times.append((t3-t2)*1000)
        out['Enumeration Method'].append('String Concatenation')
        out['Mean Processing Time (ms)'].append(round(np.mean(xs_times),8))
        out['Std Processing Time (ms)'].append(round(np.std(xs_times),8))

        hacs = [m.GetNumHeavyAtoms() for m in mols_2]

        for _ in range(2):
            #Common to string concat and reaction-based enumeration
            out['Reaction Name'].append(rxn.rxn_name)
            out['SMARTS'].append(rxn.smarts)
            out['SMARTS String Length'].append(len(rxn.smarts))
            out['Reactant 1 SMILES String Length'].append(len(rxn.example_reactant1))
            out['Mean Reactant 2 SMILES String Length'].append(round(np.mean(smiles_2_lengths),8))
            out['Std Reactant 2 SMILES String Length'].append(round(np.std(smiles_2_lengths),8))
            out['Reactant 1 HAC'].append(mols_1[0].GetNumHeavyAtoms())
            out['Mean Reactant 2 HAC'].append(round(np.mean(hacs),8))
            out['Std Reactant 2 HAC'].append(round(np.std(hacs),8))
            out['N Batches'].append(N)
            out['Batch Size'].append(BATCH_SIZE)
        print(f"""
              Reaction: {rxn.rxn_name}
              SMARTS: {rxn.smarts}
              Reactant 1 SMILES: {rxn.example_reactant1}
              Example product of string concatenation: {xs_prod_smis[0]}
              Example product of reaction: {rxn_prod_smis[0]}
              Mean string concatenation processing time: {round(np.mean(xs_times),8)}
              Mean reaction processing time: {round(np.mean(rxn_times),8)}
              Stdev string concatenation processing time: {round(np.std(xs_times),8)}
              Stdev reaction processing time:{round(np.std(rxn_times),8)}
              Mean Reactant 2 SMILES length: {round(np.mean(smiles_2_lengths),8)}
              Stdev Reactant 2 SMILES length: {round(np.std(smiles_2_lengths),8)}
              Mean Reactant 2 HAC: {round(np.mean(hacs),8)}
              Stdev Reactant 2 HAC: {round(np.std(hacs),8)}
              __________________________________________________________
              """)
    pd.DataFrame(out).to_csv('XSMILES_vs_Reaction_Speed_Benchmarking.csv',index=False)