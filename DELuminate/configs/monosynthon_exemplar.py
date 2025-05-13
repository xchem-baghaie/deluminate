COMMON_OVERWRITES = {
    'C#CCCN(CCNC(=O)OCC1c2ccccc2-c2ccccc21)CCC(=O)O':'N(C(=O)OCC1c2ccccc2-c2ccccc21)C(CC#C)C(=O)O', # Alkynyl Fmoc AA
    'C=CCCC(=O)O':'C=CCC(=O)O', # Vinyl Acid
    'CC(=O)CCCN(CCNC(=O)OC(C)(C)C)CCC(=O)O':'OC(C1CC(NC(OC(C)(C)C)=O)CC(C1)=O)=O', # Boc Amino Ketone Acid
    'CC(C)(C)OC(=O)NCCC=C[N+](=O)[O-]':'O=C(OC(C)(C)C)N1CCC=C([N+]([O-])=O)C1', # Nitroalkene
    'CC(C)(C)OC(=O)NCCC=O':'CC(C)(C)OC(=O)NCC=O', # Formyl Boc Amine
    'CC(C)(C)OC(=O)NCCN(CCNC(=O)OCC1c2ccccc2-c2ccccc21)CCC(=O)O':'CN(CC(CNC(=O)OCC1c2ccccc2-c2ccccc21)CC(=O)O)C(=O)OC(C)(C)C', # Boc Fmoc Diamino (both primary) Acid
    'CCOC(=O)CCC(Cl)=NO':'CCOC(=O)CC(Cl)=NO', # Oximinoyl Chloride Ester
    'CN(CCN(CCNC(=O)OCC1c2ccccc2-c2ccccc21)CCC(=O)O)C(=O)OC(C)(C)C':'CN(CC(CNC(=O)OCC1c2ccccc2-c2ccccc21)CC(=O)O)C(=O)OC(C)(C)C', # Boc secondary amino Fmoc primary amino Acid
    'CN(CCc1ccc(B(O)O)cc1)C(=O)OC(C)(C)C':'OB(O)C1=CCN(C(OC(C)(C)C)=O)C1', # Boc Amino Boronate
    'NCCC(=O)O':'NCC(=O)O', # Primary Amino Acid
    'NCCC1CC=CC1':'NC1CC=CC1', # Primary Amino Cyclic Alkene
    'NCCc1ccc(Br)cc1':'NCc1ccc(Br)cc1', # Primary Amino Aryl Bromide
    'O=C(O)CCC1CC=CC1':'O=C(O)C1CC=CC1', # Acid Cyclic Alkene
    'O=C(O)CCNC(=O)OCC1c2ccccc2-c2ccccc21':'O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21', # Fmoc Amino Acid
    'O=C(O)CCc1ccc(Br)cc1':'O=C(O)c1ccc(Br)cc1', # Acid Aryl Bromide
    'O=C(O)CCc1ccnc(F)c1':'O=C(O)c1ccnc(F)c1', # Acid SNAr F
    'O=CCCC(=O)O':'O=CCC(=O)O', # Formyl Acid
    'O=CCCc1ccc(F)c([N+](=O)[O-])c1':'O=Cc1ccc(F)c([N+](=O)[O-])c1', # Ortho Fluoro-Nitro Aldehyde
}
