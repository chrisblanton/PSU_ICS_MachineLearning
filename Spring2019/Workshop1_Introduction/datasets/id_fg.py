#!/usr/bin/env python
#The goal here is to id the prescence of difference functional groups
# within SMILES codes
# Smiles for 514-10-3
# CC(C)C1=CC2=CC[C@@H]3[C@@]([C@H]2CC1)(CCC[C@@]3(C)C(=O)O)C
# 
from __future__ import print_function
import os
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import FragmentCatalog
fName = os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
fparams = FragmentCatalog.FragCatParams(1,6,fName)
#print(fparams.GetNumFuncGroups())
fcat = FragmentCatalog.FragCatalog(fparams)
fcgen = FragmentCatalog.FragCatGenerator()

m = Chem.MolFromSmiles('CC(C)C1=CC2=CC[C@@H]3[C@@]([C@H]2CC1)(CCC[C@@]3(C)C(=O)O)C')
#Draw.MolToFile(m,'mymol.png')
#print(fcgen.AddFragsFromMol(m,fcat))
#print(fcat.GetEntryDescription(5))

print(list(fcat.GetEntryFuncGroupIds()))

