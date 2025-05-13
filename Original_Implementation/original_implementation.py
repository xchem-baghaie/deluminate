import pandas as pd
import numpy as np
import os
from os.path import exists
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from time import time
import random
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from PIL import Image

filterlists = False
filternumundefinedstereocenters = False
enumsubset = False
generatepropertyprofiles = True
idfmocs = False
checkwhetherstereoisomerscanphysicallyexist = False
filterforlibrarydeckinlastcycle = False
includebbswithoutvalidationdata = False

libs = [
    #'ENPL004',
    #'ENPL002',
    #'ENPL003',
    #'ENPL005',
    # 'ENPL006',
    # 'ENPL007',
    # 'ENPL008',
    'XLIB0123',
    'XLIB0124',
    'PL_2023-4',
]

valcutoff = 0.5

def main(libs):
    libinfo = pd.read_csv('Library Info.csv').fillna('')

    for lib in libs:
        print(f'__________{lib}__________\n')
        start = time()
        
        libindex = -1
        for i in range(len(libinfo)):
            if libinfo['Library'][i] == lib:
                libindex = i
                break
        if libindex == -1:
            print(f'{lib} not found in Library Info.csv')
            sys.exit()

        linkersmi = ''
        linkersmi = libinfo['Linker SMILES'][libindex]
        if linkersmi == '':
            print(f'Linker SMILES for {lib} not provided in Library Info.csv')
            sys.exit()
        
        if libinfo['Synonym'][libindex] != '':
            synonym = libinfo['Synonym'][libindex]
        else:
            synonym = lib

        dummysmiles_dct = {
            'A': {'SMILES': libinfo['Locksmiles_A'][libindex], 'Deprotection': libinfo['Deprotection_A'][libindex]},
            'B': {'SMILES': libinfo['Locksmiles_B'][libindex], 'Deprotection': libinfo['Deprotection_B'][libindex]},
            'C': {'SMILES': libinfo['Locksmiles_C'][libindex], 'Deprotection': libinfo['Deprotection_C'][libindex]},
            'D': {'SMILES': libinfo['Locksmiles_D'][libindex], 'Deprotection': libinfo['Deprotection_D'][libindex]},
            'E': {'SMILES': libinfo['Locksmiles_E'][libindex], 'Deprotection': libinfo['Deprotection_E'][libindex]},
        }

        main_dir, pre_dir, lab_dir, post_dir, individualprops_dir, propsgrid_dir = createdirs(lib)

        if filterlists:
            filterlistsfxn(linkersmi, valcutoff, dummysmiles_dct, libinfo, libindex, main_dir, pre_dir, lab_dir, post_dir)
        if filternumundefinedstereocenters:
            filternumundefinedstereocentersfxn(lab_dir, post_dir)
        if enumsubset:
            enumsubsetfxn(lib, linkersmi, main_dir, post_dir)
        
        bbfiles = os.listdir(post_dir)
        if bbfiles == []:
            print('No files were found in', post_dir)
            sys.exit()

        cycles = []

        for file in bbfiles:
            cycles.append(file[0])

        splitsizes = []

        for cycle in cycles:
            splitsizes.append(len(pd.read_csv(post_dir + cycle + '.csv')))

        libsize = 1
        for splitsize in splitsizes:
            libsize = libsize*splitsize
        
        libsizedf = pd.DataFrame({'Cycle': cycles, 'Splitsize': splitsizes, 'Library Size': [str(libsize)] + ['']*(len(cycles)-1)})
        libsizedf.to_csv(main_dir + lib + '_Split_Sizes_and_Library_Size.csv', index = False)

        if generatepropertyprofiles:
            generatepropertyprofilesfxn(lib, synonym, main_dir, individualprops_dir, propsgrid_dir)

        end = time()
        duration = int((end-start)/60)
        print(f'Processing for {lib} took {duration} minutes\n')

def createdirs(lib):
    main_dir = '/mnt/p/Discovery Chemistry/Library Synthesis/Library Planning/BB List Generation/' + lib + '/'
    pre_dir = main_dir + 'BB lists/Prefilter/'
    lab_dir = main_dir + 'BB lists/Labeled/'
    post_dir = main_dir + 'BB lists/Postfilter/'
    prof_dir = main_dir + 'Property Profiles/'
    individualprops_dir = prof_dir + 'Individual/'
    propsgrid_dir = prof_dir + 'Grid/'

    for dir in [main_dir, pre_dir, lab_dir, post_dir, prof_dir, individualprops_dir, propsgrid_dir]:
        if not exists(dir):
            os.makedirs(dir)

    return main_dir, pre_dir, lab_dir, post_dir, individualprops_dir, propsgrid_dir

def filterlistsfxn(linkersmi, valcutoff, dummysmiles_dct, libinfo, libindex, main_dir, pre_dir, lab_dir, post_dir):
    print('Filtering BB lists')
    start = time()
    bbfiles = os.listdir(pre_dir)
    if bbfiles == []:
        print('No files were found in', pre_dir)
        sys.exit()

    A_ids = []
    A_smi = []
    A_val1 = []
    A_val2 = []
    A_val3 = []
    A_val4 = []
    A_rxns = []
    A_deprs = []
    A_class = []
    B_ids = []
    B_smi = []
    B_val1 = []
    B_val2 = []
    B_val3 = []
    B_val4 = []
    B_rxns = []
    B_deprs = []
    B_class = []
    C_ids = []
    C_smi = []
    C_val1 = []
    C_val2 = []
    C_val3 = []
    C_val4 = []
    C_rxns = []
    C_deprs = []
    C_class = []
    D_ids = []
    D_smi = []
    D_val1 = []
    D_val2 = []
    D_val3 = []
    D_val4 = []
    D_rxns = []
    D_deprs = []
    D_class = []
    E_ids = []
    E_smi = []
    E_val1 = []
    E_val2 = []
    E_val3 = []
    E_val4 = []
    E_rxns = []
    E_deprs = []
    E_class = []
    A_valrxn1 = []
    B_valrxn1 = []
    C_valrxn1 = []
    D_valrxn1 = []
    E_valrxn1 = []
    A_valrxn2 = []
    B_valrxn2 = []
    C_valrxn2 = []
    D_valrxn2 = []
    E_valrxn2 = []
    A_valrxn3 = []
    B_valrxn3 = []
    C_valrxn3 = []
    D_valrxn3 = []
    E_valrxn3 = []
    A_valrxn4 = []
    B_valrxn4 = []
    C_valrxn4 = []
    D_valrxn4 = []
    E_valrxn4 = []
    A_deck = []
    B_deck = []
    C_deck = []
    D_deck = []
    E_deck = []
    
    cycles = []

    for file in bbfiles:
        print('Processing file', file)
        tmp = pd.read_csv(pre_dir + file).dropna(subset = ['SMILES']).reset_index(drop=True)
        rxncol = ''
        deprcol = ''
        for col in list(tmp.columns):
            if 'eaction' in col:
                rxncol = col
            elif 'eprotection' in col:
                deprcol = col
        if rxncol == '' or deprcol == '':
            print('Reaction or Deprotection column was not provided in BB file', pre_dir + file)
            sys.exit()
        ids = list(tmp['XBBID'])
        smi = list(tmp['SMILES'])
        rxns = list(tmp[rxncol])
        deprs = list(tmp[deprcol])
        classnames = list(tmp['CLASSNAME'])
        if 'STEP_1_RESULT' in list(tmp.columns):
            val1 = list(tmp['STEP_1_RESULT'])
            valrxn1 = list(tmp['STEP_1_ASSAY'])
        else:
            val1 = ['']*len(tmp)
            valrxn1 = ['']*len(tmp)

        if 'STEP_2_RESULT' in list(tmp.columns):
            val2 = list(tmp['STEP_2_RESULT'])
            valrxn2 = list(tmp['STEP_2_ASSAY'])
        else:
            val2 = ['']*len(tmp)
            valrxn2 = ['']*len(tmp)

        if 'STEP_3_RESULT' in list(tmp.columns):
            val3 = list(tmp['STEP_3_RESULT'])
            valrxn3 = list(tmp['STEP_3_ASSAY'])
        else:
            val3 = ['']*len(tmp)
            valrxn3 = ['']*len(tmp)

        if 'STEP_4_RESULT' in list(tmp.columns):
            val4 = list(tmp['STEP_4_RESULT'])
            valrxn4 = list(tmp['STEP_4_ASSAY'])
        else:
            val4 = ['']*len(tmp)
            valrxn4 = ['']*len(tmp)

        if 'LIBRARY_DECK_NAME' in list(tmp.columns):
            deck = list(tmp['LIBRARY_DECK_NAME'])
        else:
            deck = ['']*len(tmp)

        fmoc = Chem.MolFromSmiles('NC(OCC1C(C=CC=C2)=C2C3=C1C=CC=C3)=O')
        mols = [Chem.MolFromSmiles(x) for x in smi]
        fmocbbs = []
        fmocsmis = []
        if idfmocs:
            for i in range(len(mols)):
                print(f'Checking for Fmoc group in molecule {i+1}/{len(tmp)}', end = '\r')
                if mols[i].HasSubstructMatch(fmoc):
                    deprs[i] = 'CA DeFmoc'
                    fmocbbs.append(ids[i])
                    fmocsmis.append(smi[i])
            fmocdf = pd.DataFrame({'XBBID': fmocbbs, 'SMILES': fmocsmis})
            if len(fmocdf) > 0:
                fmocdf.to_csv(main_dir + str(len(fmocbbs)) + ' Fmoc containing BBs potentially in ' + file[:len(file)-4] + ' file.csv', index = False)
            print('')

        if file[0] == 'A':
            cycles.append('A')
            A_ids.extend(ids)
            A_smi.extend(smi)
            A_deck.extend(deck)
            A_val1.extend(val1)
            A_val2.extend(val2)
            A_val3.extend(val3)
            A_val4.extend(val4)
            A_valrxn1.extend(valrxn1)
            A_valrxn2.extend(valrxn2)
            A_valrxn3.extend(valrxn3)
            A_valrxn4.extend(valrxn4)
            A_rxns.extend(rxns)
            A_deprs.extend(deprs)
            A_class.extend(classnames)
        elif file[0] == 'B':
            cycles.append('B')
            B_ids.extend(ids)
            B_smi.extend(smi)
            B_deck.extend(deck)
            B_val1.extend(val1)
            B_val2.extend(val2)
            B_val3.extend(val3)
            B_val4.extend(val4)
            B_valrxn1.extend(valrxn1)
            B_valrxn2.extend(valrxn2)
            B_valrxn3.extend(valrxn3)
            B_valrxn4.extend(valrxn4)
            B_rxns.extend(rxns)
            B_deprs.extend(deprs)
            B_class.extend(classnames)
        elif file[0] == 'C':
            cycles.append('C')
            C_ids.extend(ids)
            C_smi.extend(smi)
            C_deck.extend(deck)
            C_val1.extend(val1)
            C_val2.extend(val2)
            C_val3.extend(val3)
            C_val4.extend(val4)
            C_valrxn1.extend(valrxn1)
            C_valrxn2.extend(valrxn2)
            C_valrxn3.extend(valrxn3)
            C_valrxn4.extend(valrxn4)
            C_rxns.extend(rxns)
            C_deprs.extend(deprs)
            C_class.extend(classnames)
        elif file[0] == 'D':
            cycles.append('D')
            D_ids.extend(ids)
            D_smi.extend(smi)
            D_deck.extend(deck)
            D_val1.extend(val1)
            D_val2.extend(val2)
            D_val3.extend(val3)
            D_val4.extend(val4)
            D_valrxn1.extend(valrxn1)
            D_valrxn2.extend(valrxn2)
            D_valrxn3.extend(valrxn3)
            D_valrxn4.extend(valrxn4)
            D_rxns.extend(rxns)
            D_deprs.extend(deprs)
            D_class.extend(classnames)
        elif file[0] == 'E':
            cycles.append('E')
            E_ids.extend(ids)
            E_smi.extend(smi)
            E_deck.extend(deck)
            E_val1.extend(val1)
            E_val2.extend(val2)
            E_val3.extend(val3)
            E_val4.extend(val4)
            E_valrxn1.extend(valrxn1)
            E_valrxn2.extend(valrxn2)
            E_valrxn3.extend(valrxn3)
            E_valrxn4.extend(valrxn4)
            E_rxns.extend(rxns)
            E_deprs.extend(deprs)
            E_class.extend(classnames)

    A_df = pd.DataFrame({'XBBID':A_ids, 'SMILES':A_smi, 'Reaction': A_rxns, 'Deprotection': A_deprs, 'Deck': A_deck, 'Classname': A_class, 'Step 1 Reaction': A_valrxn1, 'Step 1 Val Yield': A_val1, 'Step 2 Reaction': A_valrxn2, 'Step 2 Val Yield': A_val2, 'Step 3 Reaction': A_valrxn3, 'Step 3 Val Yield': A_val3, 'Step 4 Reaction': A_valrxn4, 'Step 4 Val Yield': A_val4})
    B_df = pd.DataFrame({'XBBID':B_ids, 'SMILES':B_smi, 'Reaction': B_rxns, 'Deprotection': B_deprs, 'Deck': B_deck, 'Classname': B_class, 'Step 1 Reaction': B_valrxn1, 'Step 1 Val Yield': B_val1, 'Step 2 Reaction': B_valrxn2, 'Step 2 Val Yield': B_val2, 'Step 3 Reaction': B_valrxn3, 'Step 3 Val Yield': B_val3, 'Step 4 Reaction': B_valrxn4, 'Step 4 Val Yield': B_val4})
    C_df = pd.DataFrame({'XBBID':C_ids, 'SMILES':C_smi, 'Reaction': C_rxns, 'Deprotection': C_deprs, 'Deck': C_deck, 'Classname': C_class, 'Step 1 Reaction': C_valrxn1, 'Step 1 Val Yield': C_val1, 'Step 2 Reaction': C_valrxn2, 'Step 2 Val Yield': C_val2, 'Step 3 Reaction': C_valrxn3, 'Step 3 Val Yield': C_val3, 'Step 4 Reaction': C_valrxn4, 'Step 4 Val Yield': C_val4})
    D_df = pd.DataFrame({'XBBID':D_ids, 'SMILES':D_smi, 'Reaction': D_rxns, 'Deprotection': D_deprs, 'Deck': D_deck, 'Classname': D_class, 'Step 1 Reaction': D_valrxn1, 'Step 1 Val Yield': D_val1, 'Step 2 Reaction': D_valrxn2, 'Step 2 Val Yield': D_val2, 'Step 3 Reaction': D_valrxn3, 'Step 3 Val Yield': D_val3, 'Step 4 Reaction': D_valrxn4, 'Step 4 Val Yield': D_val4})
    E_df = pd.DataFrame({'XBBID':E_ids, 'SMILES':E_smi, 'Reaction': E_rxns, 'Deprotection': E_deprs, 'Deck': E_deck, 'Classname': E_class, 'Step 1 Reaction': E_valrxn1, 'Step 1 Val Yield': E_val1, 'Step 2 Reaction': E_valrxn2, 'Step 2 Val Yield': E_val2, 'Step 3 Reaction': E_valrxn3, 'Step 3 Val Yield': E_val3, 'Step 4 Reaction': E_valrxn4, 'Step 4 Val Yield': E_val4})


    df_dct = {'A': A_df, 'B': B_df, 'C': C_df, 'D': D_df, 'E': E_df}
    rxn_dct = {}

    for cyc in df_dct.keys():
        rs = list(set(df_dct[cyc]['Reaction']))
        rs.sort()
        rxn_dct[cyc] = rs

    cycles = list(set(cycles))
    cycles.sort()
    print('\nCycles:', cycles)
    
    for cyc in rxn_dct.keys():
        if len(rxn_dct[cyc]) > 0:
            print(f'Cycle {cyc} reactions: {rxn_dct[cyc]}')

    print('')

    for cycle in cycles:
        print('######## Cycle', cycle, '########')
        df = df_dct[cycle]
        currmols = [Chem.MolFromSmiles(linkersmi)] * len(df)
        for mol in currmols:
            mol.UpdatePropertyCache()
        bbmols = [Chem.MolFromSmiles(x) for x in df['SMILES']]
        for mol in bbmols:
            mol.UpdatePropertyCache()
        df['Include'] = ['Yes']*len(df)
        df['Exclusion Rationale'] = ['']*len(df)
        df['Reactive FG'] = ['']*len(df)
        df['Incompatible FGs'] = ['']*len(df)
        df['Alert FGs'] = ['']*len(df)
        df['Warning FGs'] = ['']*len(df)


        if cycle == 'A':
            df, currmols = monosynthenumcycle(df, currmols, bbmols, valcutoff, cycle = 'A')
        elif 'A' in cycles:
            df, currmols = monosynthreactcycle(df, currmols, rxn_dct, dummysmiles_dct, cycle = 'A')
        else:
            pass

        if cycle == 'B':
            df, currmols = monosynthenumcycle(df, currmols, bbmols, valcutoff, cycle = 'B')
        elif 'B' in cycles:
            df, currmols = monosynthreactcycle(df, currmols, rxn_dct, dummysmiles_dct, cycle = 'B')
        else:
            pass

        if cycle == 'C':
            df, currmols = monosynthenumcycle(df, currmols, bbmols, valcutoff, cycle = 'C')
        elif 'C' in cycles:
            df, currmols = monosynthreactcycle(df, currmols, rxn_dct, dummysmiles_dct, cycle = 'C')
        else:
            pass

        if cycle == 'D':
            df, currmols = monosynthenumcycle(df, currmols, bbmols, valcutoff, cycle = 'D')
        elif 'D' in cycles:
            df, currmols = monosynthreactcycle(df, currmols, rxn_dct, dummysmiles_dct, cycle = 'D')
        else:
            pass

        if cycle == 'E':
            df, currmols = monosynthenumcycle(df, currmols, bbmols, valcutoff, cycle = 'E')
        elif 'E' in cycles:
            df, currmols = monosynthreactcycle(df, currmols, rxn_dct, dummysmiles_dct, cycle = 'E')
        else:
            pass
        
    
        print('Checking for duplicate monosynthon structures')

        mols = [Chem.MolFromSmiles(x) for x in df['SMILES_After_' + cycles[-1]]]
        cansmis = [Chem.MolToSmiles(x) for x in mols]
        for i in range(len(df)):
            idx = df['XBBID'][i]
            cansmi = cansmis[i]
            if cansmi != '':
                for j in range(len(df)):
                    if j > i and df['Include'][i] == 'Yes':
                        if cansmis[j] == cansmi:
                            df.loc[j, ['Include']] = 'No'
                            df.loc[j, ['Exclusion Rationale']] = f'Duplicate monosynthon enumeration of {idx}'

        print('Checking for undesirable FGs')
        
        for i in range(len(df)):
            mol = mols[i]
            excludefgs = []
            warnfgs = []
            for FG, disp in fgdisp.items():
                if disp == 'Exclude':
                    if mol.HasSubstructMatch(Chem.MolFromSmarts(fgsmarts[FG])):
                        df.loc[i, ['Include']] = 'No'
                        df.loc[i, ['Exclusion Rationale']] = 'Contains Alert Functional Group'
                        excludefgs.append(FG)
                if disp == 'Warn':
                    if mol.HasSubstructMatch(Chem.MolFromSmarts(fgsmarts[FG])):
                        warnfgs.append(FG)
            excludestring = ''
            warnstring = ''
            if len(excludefgs) > 0:
                for fg in excludefgs:
                    if excludefgs.index(fg) != len(excludefgs) - 1:
                        excludestring = excludestring + fg + '...'
                    else:
                        excludestring = excludestring + fg
            if len(warnfgs) > 0:
                for fg in warnfgs:
                    if warnfgs.index(fg) != len(warnfgs) - 1:
                        warnstring = warnstring + fg + '...'
                    else:
                        warnstring = warnstring + fg
                
            df.loc[i, ['Alert FGs']] = excludestring
            df.loc[i, ['Warning FGs']] = warnstring

        if filterforlibrarydeckinlastcycle and cycle == cycles[-1]:
            print(f'\nFiltering out BBs from cycle {cycle} according to deck assignments')
            deck = libinfo['Deck'][libindex]
            for i in range(len(df)):
                if df['Deck'][i] != deck:
                        df.loc[i, ['Include']] = 'No'
                        df.loc[i, ['Exclusion Rationale']] = f'BB is not assigned to {deck}'

        df.to_csv(lab_dir + cycle + '.csv', index = False)
        filtereddf = df[df['Include'] == 'Yes'].reset_index(drop = True)
        filtereddf.to_csv(post_dir + cycle + '.csv', index = False)

    end = time()
    enumduration = int((end-start)/60)
    print(f'BB list filtering took {enumduration} minutes\n')
    return None

def monosynthenumcycle(df, currmols, bbmols, valcutoff, cycle):
    print('Enumerating cycle', cycle)
    start = time()
    outsmi = []
    for i in range(len(df)):
        print(f'Reaction {i+1} of {len(df)}', end = '\r')
        
        if df['Classname'][i] in classnames_to_exclude:
            df.loc[i, ['Include']] = 'No'
            df.loc[i, ['Exclusion Rationale']] = 'Registered With Undesired Classname'
        
        if df['XBBID'][i].startswith('TBB'):
            df.loc[i, ['Include']] = 'No'
            df.loc[i, ['Exclusion Rationale']] = 'BBID starts with TBB'

        rxn = AllChem.ReactionFromSmarts(rxnsmarts[df['Reaction'][i]])
        if Chem.MolToSmiles(currmols[i]) != '':
            x=0
            for fg in compatible_fgs_reactant_2[df['Reaction'][i]]:
                if bbmols[i].HasSubstructMatch(Chem.MolFromSmarts(fgsmarts[fg])):
                    x+=1
                    df.loc[i, ['Reactive FG']] = fg
                if x == 0:
                    df.loc[i, ['Reactive FG']] = 'None'
                    df.loc[i, ['Include']] = 'No'
                    df.loc[i, ['Exclusion Rationale']] = 'BB does not contain a reactive functional group'
                if x > 1:
                    df.loc[i, ['Reactive FG']] = 'Multiple'
                    df.loc[i, ['Include']] = 'No'
                    df.loc[i, ['Exclusion Rationale']] = 'BB contains multiple reactive functional groups'

            x=0
            for fg in incompatible_fgs_reactant_2[df['Reaction'][i]]:
                if bbmols[i].HasSubstructMatch(Chem.MolFromSmarts(fgsmarts[fg])):
                    x+=1
                if x > 0:
                    df.loc[i, ['Incompatible FGs']] = fg
                    df.loc[i, ['Include']] = 'No'
                    df.loc[i, ['Exclusion Rationale']] = 'BB contains an incompatible functional group'

            currmols[i] = rxn.RunReactants((currmols[i], bbmols[i]))
            if len(currmols[i]) > 0:
                mol = currmols[i][0][0]
            else:
                mol = nonemol
            bbmol = bbmols[i]
            mol.UpdatePropertyCache()
            bbmol.UpdatePropertyCache()
            multireact = rxn.RunReactants((mol, bbmol))
            multireactsmi = list(set([Chem.MolToSmiles(x[0]) for x in multireact]))
            if len(multireactsmi) > 0:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'BB is capable of reacting multiple times'
            currsmi = list(set([Chem.MolToSmiles(x[0]) for x in currmols[i]]))
            if len(currsmi) != 1:
                currmols[i] = nonemol #Fail the enumeration if multiple products can form
            else:
                currmols[i] = currmols[i][0][0]
            
            currmols[i].UpdatePropertyCache()
            
            if not pd.isna(df['Deprotection'][i]):
                dep = AllChem.ReactionFromSmarts(rxnsmarts[df['Deprotection'][i]])
                currmols[i] = dep.RunReactants([currmols[i]])
                if len(currmols[i]) > 0:
                    mol = currmols[i][0][0]
                else:
                    mol = nonemol
                mol.UpdatePropertyCache()
                multireact = dep.RunReactants([mol])
                multireactsmi = list(set([Chem.MolToSmiles(x[0]) for x in multireact]))
                if len(multireactsmi) > 0:
                    df.loc[i, ['Include']] = 'No'
                    df.loc[i, ['Exclusion Rationale']] = 'BB is capable of deprotecting multiple times'
                currsmi = list(set([Chem.MolToSmiles(x[0]) for x in currmols[i]]))
                if len(currsmi) != 1:
                    currmols[i] = nonemol
                else:
                    currmols[i] = currmols[i][0][0]
            if Chem.MolToSmiles(currmols[i]) == '':
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'Failed to enumerate'
                outsmi.append('') #Failed in this step
            else:
                outsmi.append(Chem.MolToSmiles(currmols[i]))
        else:
            outsmi.append('') #Perpetuated from a previous step

        if df['Step 1 Val Yield'][i] != '':
            if df['Step 1 Val Yield'][i] < valcutoff:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'Validation failure'
        elif df['Step 1 Reaction'][i] != '':
            if not includebbswithoutvalidationdata:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'No Validation Data'                
        if df['Step 2 Val Yield'][i] != '':
            if df['Step 2 Val Yield'][i] < valcutoff**2:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'Validation failure'
        elif df['Step 2 Reaction'][i] != '':
            if not includebbswithoutvalidationdata:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'No Validation Data'   
        if df['Step 3 Val Yield'][i] != '':
            if df['Step 3 Val Yield'][i] < valcutoff**3:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'Validation failure'
        elif df['Step 3 Reaction'][i] != '':
            if not includebbswithoutvalidationdata:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'No Validation Data'   
        if df['Step 4 Val Yield'][i] != '':
            if df['Step 4 Val Yield'][i] < valcutoff**4:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'Validation failure'
        elif df['Step 4 Reaction'][i] != '':
            if not includebbswithoutvalidationdata:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'No Validation Data'   
    df['SMILES_After_' + cycle] = outsmi
    print('')
    for mol in currmols:
        mol.UpdatePropertyCache()
    end = time()
    enumduration = int((end-start)/60)
    print(f'Monosynthon enumeration for cycle {cycle} took {enumduration} minutes')
    return df, currmols

def monosynthreactcycle(df, currmols, rxn_dct, dummysmiles_dct, cycle):
    start = time()
    reactionname = rxn_dct[cycle][0] #Just use the first reaction
    if dummysmiles_dct[cycle]['SMILES'] != '':
        exemplar = Chem.MolFromSmiles(dummysmiles_dct[cycle]['SMILES'])
    else:
        exemplar = Chem.MolFromSmiles(bb_exemplars[reactionname][1])
    print('Reacting cycle', cycle, '(' + reactionname + ' with ' + Chem.MolToSmiles(exemplar) +')')
    outsmi = []
    rxn = AllChem.ReactionFromSmarts(rxnsmarts[reactionname]) 
    for i in range(len(df)):
        print(f'Reaction {i+1} of {len(df)}', end = '\r')
        if Chem.MolToSmiles(currmols[i]) != '':
            currmols[i] = rxn.RunReactants((currmols[i], exemplar))
            if len(currmols[i]) > 0:
                mol = currmols[i][0][0]
            else:
                mol = nonemol
            mol.UpdatePropertyCache()
            exemplar.UpdatePropertyCache()
            multireact = rxn.RunReactants((mol, exemplar))
            multireactsmi = list(set([Chem.MolToSmiles(x[0]) for x in multireact]))
            if len(multireactsmi) > 0:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = 'BB is capable of reacting multiple times'
            currsmi = list(set([Chem.MolToSmiles(x[0]) for x in currmols[i]]))
            if len(currsmi) != 1:
                currmols[i] = nonemol #Fail the enumeration if multiple products can form
            else:
                currmols[i] = currmols[i][0][0]
                currmols[i].UpdatePropertyCache()
                if dummysmiles_dct[cycle]['SMILES'] != '' and dummysmiles_dct[cycle]['Deprotection'] != '':
                    dep = AllChem.ReactionFromSmarts(rxnsmarts[dummysmiles_dct[cycle]['Deprotection']])
                    currmols[i] = dep.RunReactants([currmols[i]])
                    if len(currmols[i]) > 0:
                        mol = currmols[i][0][0]
                    else:
                        mol = nonemol
                    mol.UpdatePropertyCache()
                    multireact = dep.RunReactants([mol])
                    multireactsmi = list(set([Chem.MolToSmiles(x[0]) for x in multireact]))
                    if len(multireactsmi) > 0:
                        df.loc[i, ['Include']] = 'No'
                        df.loc[i, ['Exclusion Rationale']] = 'BB is capable of deprotecting multiple times'
                    currsmi = list(set([Chem.MolToSmiles(x[0]) for x in currmols[i]]))
                    if len(currsmi) != 1:
                        currmols[i] = nonemol #Fail the enumeration if multiple products can form
                    else:
                        currmols[i] = currmols[i][0][0]
            if Chem.MolToSmiles(currmols[i]) == '':
                outsmi.append('') #Failed in this step
            else:
                outsmi.append(Chem.MolToSmiles(currmols[i]))
        else:
            outsmi.append('') #Perpetuated from a previous step
        if Chem.MolToSmiles(currmols[i]) == '':
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Failed to Enumerate']] = 'Yes'
    df['SMILES_After_' + cycle] = outsmi
    print('')
    for mol in currmols:
        mol.UpdatePropertyCache()
    end = time()
    enumduration = int((end-start)/60)
    print(f'Monosynthon reaction for cycle {cycle} took {enumduration} minutes')
    return df, currmols

def filternumundefinedstereocentersfxn(lab_dir, post_dir):
    print('Filtering number of undefined stereocenters')
    start = time()
    bbfiles = os.listdir(lab_dir)
    if bbfiles == []:
        print('No files were found in', lab_dir)
        sys.exit()

    cycles = []

    for file in bbfiles:
        cycles.append(file[0])

    for cycle in cycles:
        print('Calculating number of undefined stereocenters for cycle', cycle)
        df = pd.read_csv(lab_dir + cycle + '.csv')
        df.fillna('', inplace=True)
        mols = [Chem.MolFromSmiles(x) for x in list(df['SMILES_After_' + cycles[-1]])]
        df['Num Possible Stereoisomers'] = [0]*len(df)
        for i in range(len(df)):
            print(f'Molecule {i+1}/{len(df)}', end = '\r')
            if Chem.MolToSmiles(mols[i]) != '':
                undefinedisom = [Chem.MolToSmiles(x) for x in tuple(EnumerateStereoisomers(mols[i], 
                options=StereoEnumerationOptions(tryEmbedding=checkwhetherstereoisomerscanphysicallyexist, onlyUnassigned=True)))]
                if len(undefinedisom) > 0:
                    numpossstereoisom = len(undefinedisom)
                else:
                    numpossstereoisom = 1
                df.loc[i, 'Num Possible Stereoisomers'] = numpossstereoisom
                if numpossstereoisom > 4:
                    df.loc[i, ['Include']] = 'No'
                    df.loc[i, ['Exclusion Rationale']] = 'Mixture of >4 stereoisomers'
        df.to_csv(lab_dir + cycle + '.csv', index = False)
        filtereddf = df[df['Include'] == 'Yes'].reset_index(drop = True)
        filtereddf.to_csv(post_dir + cycle + '.csv', index = False)

    if len(cycles) == 2:
        opts = [(4,1), (1,4), (2,2)]
    elif len(cycles) == 3:
        opts = [(4,1,1), (1,4,1), (1,1,4), (2,2,1), (2,1,2), (1,2,2)]
    elif len(cycles) == 4:
        opts = [(4,1,1,1), (1,4,1,1), (1,1,4,1), (1,1,1,4), (2,2,1,1), (2,1,2,1), (2,1,1,2), (1,2,2,1), (1,2,1,2), (1,1,2,2)]

    undefinedcountsdct = {}
    for cycle in cycles:
        ls = list(pd.read_csv(post_dir + cycle + '.csv')['Num Possible Stereoisomers'])
        dct = {}
        for ct in list(set(ls)):
            dct[ct] = ls.count(ct)
        undefinedcountsdct[cycle] = dct
    
    print('Undefined Counts dictionary:', undefinedcountsdct)

    optresults = {}
    
    for opt in opts:
        filteredpct = 1
        for i in range(len(cycles)):
            bbct = 0
            allbbct = 0
            for key in undefinedcountsdct[cycles[i]].keys(): #key = undefined stereocenter counts
                if key <= 4:
                    allbbct += undefinedcountsdct[cycles[i]][key]
                if key <= opt[i]: 
                    bbct += undefinedcountsdct[cycles[i]][key]
            filteredpct = filteredpct * (bbct/allbbct)
        optresults[opt] = filteredpct
    
    print(optresults)

    best = max(list(optresults.values()))

    for k, v in optresults.items():
        if v == best:
            bestopt = k
    
    print('Best Stereochemistry Distribution Option:', bestopt)

    for cycle in cycles:
        maxallowedstereocenters = bestopt[cycles.index(cycle)]
        df = pd.read_csv(lab_dir + cycle + '.csv')
        df.fillna('', inplace=True)
        for i in range(len(df)):
            if df['Num Possible Stereoisomers'][i] > maxallowedstereocenters:
                df.loc[i, ['Include']] = 'No'
                df.loc[i, ['Exclusion Rationale']] = f'Too many possible stereoisomers ({str(maxallowedstereocenters)} allowed for cycle {cycle})'
        df.to_csv(lab_dir + cycle + '.csv', index = False)
        filtereddf = df[df['Include'] == 'Yes'].reset_index(drop = True)
        filtereddf.to_csv(post_dir + cycle + '.csv', index = False)
    end = time()
    stereochemduration = int((end-start)/60)
    print(f'Stereochemistry filtering took {stereochemduration} minutes\n')
    return None

def enumsubsetfxn(lib, linkersmi, main_dir, post_dir):
    print('Enumerating Instances')
    start = time()
    bbfiles = os.listdir(post_dir)
    if bbfiles == []:
        print('No files were found in', post_dir)
        sys.exit()

    cycles = []

    for file in bbfiles:
        cycles.append(file[0])
    
    cycles = list(set(cycles))
    cycles.sort()

    splitsizes = []
    collist = ['Library']
    for cycle in cycles:
        splitsizes.append(len(pd.read_csv(post_dir + cycle + '.csv')))
        collist.append('BB' + cycle + ' ID')
    
    for cycle in cycles:
        collist.append('BB' + cycle + ' SMILES')

    for cycle in cycles:
        collist.append('BB' + cycle + ' Reaction')
        collist.append('BB' + cycle + ' Deprotection')

    for cycle in cycles:
        collist.append('BB' + cycle + ' Substructure Warnings')

    collist.append('Enumerated SMILES')

    df = pd.DataFrame(columns = collist)

    for splitsize in splitsizes:
        if splitsize < 1:
            emptycycle = cycles[splitsizes.index(splitsize)]
            print(f'Cycle {emptycycle} is empty. Exiting...')
            sys.exit()
    maxsplitsize = max(splitsizes)
    numinstances = max([maxsplitsize, 3000])
    print(f'{numinstances} instances will be enumerated')

    df['Library'] = [lib] * numinstances
    print('Cycles:', cycles)
    for cycle in cycles:
        k = 0
        cycdf = pd.read_csv(post_dir + cycle + '.csv')
        inds = list(cycdf.index.values)
        shufids = []
        shufsmi = []
        shufrxn = []
        shufdep = []
        shufwarn = []
        while k < numinstances:
            random.seed(k) #Different seed every iteration, but overall enumeration is reproducible
            shufinds = random.sample(inds, len(inds))
            for ind in shufinds:
                if k >= numinstances:
                    break
                else:
                    k += 1
                    shufids.append(cycdf['XBBID'][ind])
                    shufsmi.append(cycdf['SMILES'][ind])
                    shufrxn.append(cycdf['Reaction'][ind])
                    shufdep.append(cycdf['Deprotection'][ind])
                    shufwarn.append(cycdf['Warning FGs'][ind])
        df['BB' + cycle + ' ID'] = shufids
        df['BB' + cycle + ' SMILES'] = shufsmi
        df['BB' + cycle + ' Reaction'] = shufrxn
        df['BB' + cycle + ' Deprotection'] = shufdep
        df['BB' + cycle + ' Substructure Warnings'] = shufwarn

    currmols = [Chem.MolFromSmiles(linkersmi)] * numinstances

    for cycle in cycles:
        if cycle == 'A':
            bbmols = [Chem.MolFromSmiles(x) for x in df['BBA SMILES']]
            currmols = instanceenumeratecycle(df, currmols, bbmols, cycle = 'A')
        else:
            pass

        if cycle == 'B':
            bbmols = [Chem.MolFromSmiles(x) for x in df['BBB SMILES']]
            currmols = instanceenumeratecycle(df, currmols, bbmols, cycle = 'B')
        else:
            pass

        if cycle == 'C':
            bbmols = [Chem.MolFromSmiles(x) for x in df['BBC SMILES']]
            currmols = instanceenumeratecycle(df, currmols, bbmols, cycle = 'C')
        else:
            pass

        if cycle == 'D':
            bbmols = [Chem.MolFromSmiles(x) for x in df['BBD SMILES']]
            currmols = instanceenumeratecycle(df, currmols, bbmols, cycle = 'D')
        else:
            pass

        if cycle == 'E':
            bbmols = [Chem.MolFromSmiles(x) for x in df['BBE SMILES']]
            currmols = instanceenumeratecycle(df, currmols, bbmols, cycle = 'E')
        else:
            pass
    
    numfailures = 0
    for i in range(len(df)):
        if df['Enumerated SMILES'][i] == '':
            numfailures += 1

    if numfailures == 0:
        df.to_csv(main_dir + lib + '_' + str(numinstances) + '_Enumerated_Instances.csv', index = False)
        df.head(3000).to_csv(main_dir + lib + '_3k_Enumerated_Instances_Properties.csv', index = False)
    else:
        df.to_csv(main_dir + str(numfailures) + '_FAILURES_' + lib + '_' + str(numinstances) + '_Enumerated_Instances.csv', index = False)
        df.head(3000).to_csv(main_dir + str(numfailures) + '_FAILURES_' + lib + '_3k_Enumerated_Instances_Properties.csv', index = False)
    end = time()
    duration = int((end-start)/60)
    print(f'Instance enumeration took {duration} minutes\n')
    return None

def instanceenumeratecycle(df, currmols, bbmols, cycle):
    print('Enumerating instances for cycle', cycle)
    outsmi = []
    for i in range(len(df)):
        print(f'Reaction {i+1} of {len(df)}', end = '\r')
        rxn = AllChem.ReactionFromSmarts(rxnsmarts[df['BB' + cycle + ' Reaction'][i]])
        if Chem.MolToSmiles(currmols[i]) != '':
            currmols[i] = rxn.RunReactants((currmols[i], bbmols[i]))
            currsmi = list(set([Chem.MolToSmiles(x[0]) for x in currmols[i]]))
            if len(currsmi) != 1:
                currmols[i] = nonemol #Fail the enumeration if multiple products can form
            else:
                currmols[i] = currmols[i][0][0]
            
            currmols[i].UpdatePropertyCache()
            
            if not pd.isna(df['BB' + cycle + ' Deprotection'][i]):
                dep = AllChem.ReactionFromSmarts(rxnsmarts[df['BB' + cycle + ' Deprotection'][i]])
                currmols[i] = dep.RunReactants([currmols[i]])
                currsmi = list(set([Chem.MolToSmiles(x[0]) for x in currmols[i]]))
                if len(currsmi) != 1:
                    currmols[i] = nonemol
                else:
                    currmols[i] = currmols[i][0][0]
        else:
            print('Enumeration failed in a previous cycle!')
        
        outsmi.append(Chem.MolToSmiles(currmols[i]))
    df['Enumerated SMILES'] = outsmi
    print('')
    for mol in currmols: 
        if Chem.MolToSmiles(mol) != '':
            mol.UpdatePropertyCache()
    return currmols

def generatepropertyprofilesfxn(lib, synonym, main_dir, individualprops_dir, propsgrid_dir):
    print('Generating property profiles\n')
    start = time()
    filenames = os.listdir(main_dir)
    enumfile = ''
    x = 0
    for filename in filenames:
        if 'roperties.csv' in filename:
            x += 1
            enumfile = main_dir + filename
    if enumfile == '':
        print('Enumeration file not found. Exiting...')
        sys.exit()
    if x > 1:
        print('Multiple files with filenames ending in "Properties.csv" were found. Exiting...')
        sys.exit()

    df = pd.read_csv(enumfile, low_memory = False)
    for prop in props.keys():
        fig, ax = plt.subplots()
        data = df[props[prop]].tolist()
        if prop in ['pKa', 'pKb']:
            data = [float(x) for x in data if x != 'No value calculated' and x != 'Calculation not licensed']
        bins = bin_dct[prop]
        if prop != 'Lipinski Compliance':
            ax = sns.histplot(data, bins=bins, kde = False, color = color_dct[prop])
            ax.set_xticks(bin_dct[prop])
            plt.title(synonym + ' ' + prop + ' Distribution', fontsize=18)
            if prop in ['pKa', 'pKb']:
                plt.title(synonym + ' ' + prop + ' Distribution (Ionizable only)', fontsize=15)
            plt.yticks([])
            ax.get_yaxis().set_visible(False)
            fig.savefig(individualprops_dir + '/' + synonym + '_' + prop + '.png', bbox_inches = 'tight')
        else:
            unique_data = sorted(list(set(data)), reverse=True)
            ls = ['PASS', 'EXCEPTIONS: 1', 'EXCEPTIONS: 2', 'EXCEPTIONS: 3', 'EXCEPTIONS: 4']
            unique_data = sorted(list(set(unique_data).intersection(ls)))
            del unique_data[-1]
            unique_data = ['PASS'] + unique_data
            percentages = {}
            for outcome in unique_data:
                percentages[outcome] = data.count(outcome)/len(data)
            labels = percentages.keys()
            pc = percentages.values()
            colors = ['#2ecc71', '#ffc3a0', '#f08080', '#808080', '#000000']
            plt.pie(pc, labels=labels, colors = colors)#, autopct='%.0f%%')
            plt.title(synonym + ' Lipinski Ro5 Compliance', fontsize=18)
            fig.savefig(individualprops_dir  + synonym + '_' + prop + '.png', bbox_inches = 'tight')  
    mw = Image.open(individualprops_dir + synonym + '_MW.png')
    clogp = Image.open(individualprops_dir + synonym + '_cLogP.png')
    clogd = Image.open(individualprops_dir + synonym + '_cLogD (pH 7.4).png')
    fsp3 = Image.open(individualprops_dir + synonym + '_Fsp3.png')
    hbd = Image.open(individualprops_dir + synonym + '_HBD.png')
    hba = Image.open(individualprops_dir + synonym + '_HBA.png')
    tpsa = Image.open(individualprops_dir + synonym + '_tPSA.png')
    rotb = Image.open(individualprops_dir + synonym + '_RotB.png')
    x = mw.size[0]
    y = mw.size[1]
    grid = Image.new('RGB',(4*x, 2*y), (250,250,250))
    grid.paste(mw, (0, 0))
    grid.paste(clogp, (x, 0))
    grid.paste(clogd, (2*x, 0))
    grid.paste(fsp3, (3*x, 0))
    grid.paste(hbd, (0, y))
    grid.paste(hba, (x, y))
    grid.paste(tpsa, (2*x, y))
    grid.paste(rotb, (3*x, y))
    grid.save(propsgrid_dir + lib + '_grid.png')
    plt.close('all')
    return None

props = {
    'MW': 'CA_MW',
    'cLogP': 'CA_cLogP',
    'cLogD (pH 7.4)': 'CA_cLogD_7.4',
    'Fsp3': 'CA_FSP3',
    'HBD': 'CA_HBD',
    'HBA': 'CA_HBA',
    'tPSA': 'CA_tPSA',
    'RotB': 'CA_ROT',
    'pKa': 'CA_PKA',
    'pKb': 'CA_PKB',
    'Ring Count': 'CA_RING_COUNT',
    'Aromatic Ring Count': 'CA_AROMATIC_RING_COUNT',
    'Heavy Atom Count': 'CA_HAC',
    'Lipinski Compliance': 'CA_RULEOF5',
}

bin_dct = {
    'MW': range(250, 800, 50),
    'cLogP': range(-5, 10, 2),
    'cLogD (pH 7.4)': range(-5, 10, 2),
    'Fsp3': [x for x in np.arange(0, 1.1, 0.1)],
    'HBD': range(0, 8),
    'HBA': range(0, 16, 2),
    'tPSA': range(0, 220, 20),
    'RotB': range(0, 22, 2),
    'pKa': [x for x in np.arange(-6., 21., 3.)],
    'pKb': [x for x in np.arange(-6., 21., 3.)],
    'Ring Count': range(0, 11),
    'Aromatic Ring Count': range(0, 11),
    'Heavy Atom Count': range(10, 90, 10),
    'Lipinski Compliance': None,
}

color_dct = {
    'MW': '#CC99FF',
    'cLogP': '#5C89E4',
    'cLogD (pH 7.4)': '#AFD5AF',
    'Fsp3': '#9E389E',
    'HBD': '#03A32B',
    'HBA': '#6FDC6F',
    'tPSA': '#265F97',
    'RotB': '#67B1F9',
    'pKa': None,
    'pKb': None,
    'Ring Count': None,
    'Aromatic Ring Count': None,
    'Heavy Atom Count': None,
    'Lipinski Compliance': None,
}

patch_dct = {
    'MW': 5,
    'cLogP': 5,
    'cLogD (pH 7.4)': 5,
    'Fsp3': 1,
    'HBD': 5,
    'HBA': 5,
    'tPSA': 7,
    'RotB': 5,
    'pKa': False,
    'pKb': False,
    'Ring Count': False,
    'Aromatic Ring Count': False,
    'Heavy Atom Count': False,
    'Lipinski Compliance': False,
}

classnames_to_exclude = [
    'AZ_DEL Core',
    'Biotin',
    'Linkers',
    'Other',
    'Covalent Warheads',
    'TechDev_BBs',
]

fg_df = pd.read_csv('Functional Groups.csv')
rxnsmarts_df = pd.read_csv('Reaction SMARTS.csv')

fgsmarts = {}
fgdisp = {}
for i in range(len(fg_df)):
    fgsmarts[fg_df['Description'][i]] = fg_df['SMARTS'][i]
    fgdisp[fg_df['Description'][i]] = fg_df['Library BB List Filtering Rule Set'][i]

rxnsmarts = {}
for i in range(len(rxnsmarts_df)):
    rxnsmarts[rxnsmarts_df['Reaction'][i]] = rxnsmarts_df['SMARTS'][i]

rxns = {
    '1,2,4-Oxadiazole Formation': 'CA 124 Oxadiazole',
    '1,3,4-Oxadiazole Formation': 'CA 134 Oxadiazole',
    'Acylation': 'CA Acylation',
    'Alkene Oxidation and Ring Expansion': 'CA Primary Amine Oxidative Ring Expansion',
    'Boc Deprotection': 'CA DeBoc',
    'Ester Hydrolysis': 'CA Ester Hydrolysis',
    'Fmoc Deprotection': 'CA DeFmoc',
    'Reductive Alkylation': 'CA Reductive Alkylation',
    'Reductive Amination': 'CA Reductive Amination',
    'Reverse Amidation': 'CA Reverse Amidation',
    'SNAr (Heteroaryl Halide BBs)': 'CA SNAr with Electrophiles',
    'SNAr (Nucleophilic BBs)': 'CA SNAr with Nucleophiles',
    'Sulfonamidation': 'CA Sulfonamidation',
    'Suzuki-sSPhos': 'CA Suzuki',
    'Urea Formation (Isocyanates)': 'CA Ureation Isocyanates',
    'CA 124 Oxadiazole': 'CA 124 Oxadiazole',
    'CA 134 Oxadiazole': 'CA 134 Oxadiazole',
    'CA Acylation': 'CA Acylation',
    'CA Primary Amine Oxidative Ring Expansion': 'CA Primary Amine Oxidative Ring Expansion',
    'CA DeBoc': 'CA DeBoc',
    'CA Ester Hydrolysis': 'CA Ester Hydrolysis',
    'CA DeFmoc': 'CA DeFmoc',
    'CA Reductive Alkylation': 'CA Reductive Alkylation',
    'CA Reductive Amination': 'CA Reductive Amination',
    'CA Reverse Amidation': 'CA Reverse Amidation',
    'CA SNAr with Electrophiles': 'CA SNAr with Electrophiles',
    'CA SNAr with Nucleophiles': 'CA SNAr with Nucleophiles',
    'CA Sulfonamidation': 'CA Sulfonamidation',
    'CA Suzuki': 'CA Suzuki',
    'CA Ureation Isocyanates': 'CA Ureation Isocyanates',
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': 'ChemAxon-UreaFormation-SecAmine-PrimAmine',
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': 'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182',
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': 'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182',
    'ChemAxon-ArylBr-to-Benzylamine': 'ChemAxon-ArylBr-to-Benzylamine',
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW': 'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW',
    'ChemAxon-PiperidoneFormation-XPL0200': 'ChemAxon-PiperidoneFormation-XPL0200',
    'ChemAxon-Acylation-Isoxazole-XPL0182':'ChemAxon-Acylation-Isoxazole-XPL0182',
}

num_reactant_dict = {
    'CA 124 Oxadiazole': 2,
    'CA 134 Oxadiazole': 2,
    'CA Acylation': 2,
    'CA Primary Amine Oxidative Ring Expansion': 2,
    'CA DeBoc': 1,
    'CA Ester Hydrolysis': 1,
    'CA DeFmoc': 1,
    'CA Reductive Alkylation': 2,
    'CA Reductive Amination': 2,
    'CA Reverse Amidation': 2,
    'CA SNAr with Electrophiles': 2,
    'CA SNAr with Nucleophiles': 2,
    'CA Sulfonamidation': 2,
    'CA Suzuki': 2,
    'CA Ureation Isocyanates': 2,
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': 2,
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': 2,
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': 2,
    'ChemAxon-ArylBr-to-Benzylamine': 1,
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW':2,
    'ChemAxon-PiperidoneFormation-XPL0200': 1,   
    'ChemAxon-Acylation-Isoxazole-XPL0182': 1, 
}

compatible_fgs_reactant_1 = {
    'CA 124 Oxadiazole': ['Secondary Aliphatic or Aromatic Amine'],
    'CA 134 Oxadiazole': ['Secondary Aliphatic or Aromatic Amine'],
    'CA Acylation': ['Primary or secondary aliphatic or aromatic amine'],
    'CA Primary Amine Oxidative Ring Expansion': ['Alkene (Cyclic)'],
    'CA DeBoc': ['Boc Protected Amine'],
    'CA Ester Hydrolysis': ['Ester'],
    'CA DeFmoc': ['Fmoc Protected Amine'],
    'CA Reductive Alkylation': ['Secondary Aliphatic Amine'],
    'CA Reductive Amination': ['Aldehyde'],
    'CA Reverse Amidation': ['Carboxylic Acid'],
    'CA SNAr with Electrophiles': ['Primary or secondary aliphatic amine'],
    'CA SNAr with Nucleophiles': ['SNAr F or Cl'],
    'CA Sulfonamidation': ['Primary or secondary aliphatic or aromatic amine'],
    'CA Suzuki': ['Aryl Cl Br or I'],
    'CA Ureation Isocyanates': ['Primary or secondary aliphatic or aromatic amine'],
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': ['Secondary Aliphatic Amine'],
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': ['Isoxazole'],
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': ['Isoxazole'],
    'ChemAxon-ArylBr-to-Benzylamine': ['Aryl Bromide'],
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW': ['Ketone', 'Aldehyde'],
    'ChemAxon-PiperidoneFormation-XPL0200': ['Fmoc Protected Primary Amine'],
    'ChemAxon-Acylation-Isoxazole-XPL0182': ['Primary or secondary aliphatic or aromatic amine'],
}

incompatible_fgs_reactant_1 = {
    'CA 124 Oxadiazole': ['Poly-Amine', 'Nitrile', 'Hydrazide', 'Hydrazine'],
    'CA 134 Oxadiazole': ['Poly-Amine', 'Hydrazide', 'Hydrazine'],
    'CA Acylation': ['Poly-Amine', 'Hydrazide', 'Hydrazine', 'Amidine'],
    'CA Primary Amine Oxidative Ring Expansion': ['Alkene (Acyclic)', 'Alkyne', 'Poly-Alkene', 'Allene', 'Primary secondary or tertiary aliphatic or aromatic amine'],
    'CA DeBoc': ['Poly-Boc Amine', 'Fmoc Protected Amine'],
    'CA Ester Hydrolysis': ['Poly-Ester', 'Fmoc Protected Amine'],
    'CA DeFmoc': ['Poly-Fmoc Amine'],
    'CA Reductive Alkylation': ['Poly-Amine', 'Hydrazide', 'Hydrazine', 'Primary Aliphatic or Aromatic Amine', 'Secondary Aromatic Amine'],
    'CA Reductive Amination': ['Poly-Aldehyde'],
    'CA Reverse Amidation': ['Poly-Carboxylic Acid', 'Primary or secondary aliphatic amine', 'Primary Aromatic Amine', 'Alkyl halide', 'Acyl Halide'],
    'CA SNAr with Electrophiles': ['Poly-Amine', 'Hydrazide', 'Hydrazine'],
    'CA SNAr with Nucleophiles': ['Poly-SNAr F or Cl', 'Alkyl halide', 'Acyl Halide'],
    'CA Sulfonamidation': ['Poly-Amine', 'Hydrazide', 'Hydrazine'],
    'CA Suzuki': ['Poly-Aryl Cl Br or I'],
    'CA Ureation Isocyanates': ['Poly-Amine'],
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': ['Poly-Amine'],
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': ['Poly-Isoxazole'],
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': ['Poly-Isoxazole'],
    'ChemAxon-ArylBr-to-Benzylamine': ['Boc Protected Amine', 'Aryl Chloride', 'Aryl Iodide'],
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW': ['Poly-Ketone', 'Poly-Aldehyde'],
    'ChemAxon-PiperidoneFormation-XPL0200': ['Primary or secondary aliphatic or aromatic amine', 'Poly-Fmoc Amine', 'Hydrazine', 'Hydrazide'],
    'ChemAxon-Acylation-Isoxazole-XPL0182': ['Poly-Amine', 'Hydrazide', 'Hydrazine', 'Amidine'],
}

compatible_fgs_reactant_2 = {
    'CA 124 Oxadiazole': ['Carboxylic Acid'],
    'CA 134 Oxadiazole': ['Aldehyde'],
    'CA Acylation': ['Carboxylic Acid', 'Acyl Halide', 'NHS Ester'],
    'CA Primary Amine Oxidative Ring Expansion': ['Primary or secondary aliphatic amine'],
    'CA DeBoc': [],
    'CA Ester Hydrolysis': [],
    'CA DeFmoc': [],
    'CA Reductive Alkylation': ['Aldehyde'],
    'CA Reductive Amination': ['Primary or secondary aliphatic or aromatic amine'],
    'CA Reverse Amidation': ['Carboxylic Acid'],
    'CA SNAr with Electrophiles': ['SNAr F or Cl'],
    'CA SNAr with Nucleophiles': ['Primary or secondary aliphatic amine or phenol'],
    'CA Sulfonamidation': ['Sulfonyl Chloride'],
    'CA Suzuki': ['Boronate (acid or ester)'],
    'CA Ureation Isocyanates': ['Isocyanate'],
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': ['Primary Aliphatic Amine'],
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': ['Hydrazine'],
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': ['Primary Aliphatic Amine'],
    'ChemAxon-ArylBr-to-Benzylamine': [],
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW': ['Primary or secondary aliphatic or aromatic amine'],
    'ChemAxon-PiperidoneFormation-XPL0200': [],
    'ChemAxon-Acylation-Isoxazole-XPL0182': [],
}

incompatible_fgs_reactant_2 = { 
    'CA 124 Oxadiazole': ['Poly-Carboxylic Acid'],
    'CA 134 Oxadiazole': ['Poly-Aldehyde'],
    'CA Acylation': ['Poly-Carboxylic Acid'],
    'CA Primary Amine Oxidative Ring Expansion': ['Poly-Amine', 'Hydrazine', 'Hydrazide'],
    'CA DeBoc': [],
    'CA Ester Hydrolysis': [],
    'CA DeFmoc': [],
    'CA Reductive Alkylation': ['Poly-Aldehyde'],
    'CA Reductive Amination': ['Poly-Amine', 'Hydrazine', 'Hydrazide'],
    'CA Reverse Amidation': ['Poly-Amine', 'Hydrazine', 'Hydrazide'],
    'CA SNAr with Electrophiles': ['Poly-SNAr F or Cl'],
    'CA SNAr with Nucleophiles': ['Poly-Amine', 'Poly-Phenol', 'Primary or secondary aliphatic amine and phenol'],
    'CA Sulfonamidation': ['Poly-Sulfonyl Chloride'],
    'CA Suzuki': ['Poly-Boronate (acid or ester)'],
    'CA Ureation Isocyanates': ['Poly-Isocyanate'],
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': ['Poly-Amine', 'Hydrazine', 'Hydrazide'],
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': ['Hydrazide'],
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': ['Hydrazine', 'Hydrazide', 'Secondary Aliphatic Amine', 'Secondary Aromatic Amine'],
    'ChemAxon-ArylBr-to-Benzylamine': [],
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW': ['Poly-Amine', 'Hydrazine', 'Hydrazide'],
    'ChemAxon-PiperidoneFormation-XPL0200': [],
    'ChemAxon-Acylation-Isoxazole-XPL0182': [],
}

bb_exemplars = { 
    'CA 124 Oxadiazole': ('CN', 'CC(=O)O'),
    'CA 134 Oxadiazole': ('CN', 'CC(=O)'),
    'CA Acylation': ('CN', 'CC(=O)O'),
    'CA Primary Amine Oxidative Ring Expansion': ('C1=CCCC1', 'CN'),
    'CA DeBoc': ('CNC(OC(C)(C)C)=O'),
    'CA Ester Hydrolysis': ('COC(=O)C'),
    'CA DeFmoc': ('CNC(OCC1C(C=CC=C2)=C2C3=C1C=CC=C3)=O'),
    'CA Reductive Alkylation': ('CNC', 'CC(=O)'),
    'CA Reductive Amination': ('CC=O', 'CN'),
    'CA Reverse Amidation': ('CC(=O)O', 'CN'),
    'CA SNAr with Electrophiles': ('CN', 'FC1=NC=CC=C1'),
    'CA SNAr with Nucleophiles': ('FC1=NC=CC=C1', 'CN'),
    'CA Sulfonamidation': ('CN', 'CS(=O)(=O)Cl'),
    'CA Suzuki': ('c1ccccc1Br', 'B(O)(O)c1ccccc1'),
    'CA Ureation Isocyanates': ('CN', 'CN=C=O'),
    'ChemAxon-UreaFormation-SecAmine-PrimAmine': ('CNC', 'CN'),
    'ChemAxon-PyrazoleFormation-Hydrazines-XPL0182': ('CN(C)C(C1=NOC=C1)=O', 'CNN'),
    'ChemAxon-PyrazoleFormation-PrimaryAmines-XPL0182': ('CN(C)C(C1=NOC=C1)=O', 'CN'),
    'ChemAxon-ArylBr-to-Benzylamine': ('c1ccccc1Br'),
    'ChemAxon-ReductiveAmination-KetonesAndAldehydes-RW':('CC(=O)C', 'CN'),
    'ChemAxon-PiperidoneFormation-XPL0200': ('CN'),
    'ChemAxon-Acylation-Isoxazole-XPL0182': ('CN'), 
}

for rxn in incompatible_fgs_reactant_2.keys():
    for x in compatible_fgs_reactant_1[rxn]:
        if x not in incompatible_fgs_reactant_2[rxn] and num_reactant_dict[rxn] == 2:
            incompatible_fgs_reactant_2[rxn].append(x)

nonemol = Chem.MolFromSmiles('')

if __name__ == '__main__':
    main(libs)