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
#From contrib ifg.py, there must be another way to use this
from collections import namedtuple

def merge(mol, marked, aset):
    bset = set()
    for idx in aset:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            jdx = nbr.GetIdx()
            if jdx in marked:
                marked.remove(jdx)
                bset.add(jdx)
    if not bset:
        return
    merge(mol, marked, bset)
    aset.update(bset)

# atoms connected by non-aromatic double or triple bond to any heteroatom
# c=O should not match (see fig1, box 15).  I think using A instead of * should sort that out?
PATT_DOUBLE_TRIPLE = Chem.MolFromSmarts('A=,#[!#6]')
# atoms in non aromatic carbon-carbon double or triple bonds
PATT_CC_DOUBLE_TRIPLE = Chem.MolFromSmarts('C=,#C')
# acetal carbons, i.e. sp3 carbons connected to tow or more oxygens, nitrogens or sulfurs; these O, N or S atoms must have only single bonds
PATT_ACETAL = Chem.MolFromSmarts('[CX4](-[O,N,S])-[O,N,S]')
# all atoms in oxirane, aziridine and thiirane rings
PATT_OXIRANE_ETC = Chem.MolFromSmarts('[O,N,S]1CC1')

PATT_TUPLE = (PATT_DOUBLE_TRIPLE, PATT_CC_DOUBLE_TRIPLE, PATT_ACETAL, PATT_OXIRANE_ETC)

def identify_functional_groups(mol):
    marked = set()
#mark all heteroatoms in a molecule, including halogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6,1): # would we ever have hydrogen?
            marked.add(atom.GetIdx())

#mark the four specific types of carbon atom
    for patt in PATT_TUPLE:
        for path in mol.GetSubstructMatches(patt):
            for atomindex in path:
                marked.add(atomindex)

#merge all connected marked atoms to a single FG
    groups = []
    while marked:
        grp = set([marked.pop()])
        merge(mol, marked, grp)
        groups.append(grp)

#extract also connected unmarked carbon atoms
    ifg = namedtuple('IFG', ['atomIds', 'atoms', 'type'])
    ifgs = []
    for g in groups:
        uca = set()
        for atomidx in g:
            for n in mol.GetAtomWithIdx(atomidx).GetNeighbors():
                if n.GetAtomicNum() == 6:
                    uca.add(n.GetIdx())
        ifgs.append(ifg(atomIds=tuple(list(g)), atoms=Chem.MolFragmentToSmiles(mol, g, canonical=True), type=Chem.MolFragmentToSmiles(mol, g.union(uca),canonical=True)))
    return ifgs
# End of ifg.py from contrib

#print(RDConfig.RDDataDir)
fName = os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
fparams = FragmentCatalog.FragCatParams(1,6,fName)
#print(fparams.GetNumFuncGroups())
fcat = FragmentCatalog.FragCatalog(fparams)
fcgen = FragmentCatalog.FragCatGenerator()



# 514-10-3
#m = Chem.MolFromSmiles('CC(C)C1=CC2=CC[C@@H]3[C@@]([C@H]2CC1)(CCC[C@@]3(C)C(=O)O)C')
m = Chem.MolFromSmiles('CC(=O)O')
# trans-Cinnamic Acid
#m = Chem.MolFromSmiles('c1ccc(cc1)/C=C/C(=O)O')
# ethanol
#m = Chem.MolFromSmiles('CCO')
# diethyl ether
#m = Chem.MolFromSmiles('CCOCC')
#Draw.MolToFile(m,'mymol.png')
#print(fcgen.AddFragsFromMol(m,fcat))
#print(fcat.GetEntryDescription(5))

#print(fcat.GetNumEntries())
#for i in range(fcat.GetNumEntries()):
#    funcGroupList = list(fcat.GetEntryFuncGroupIds(i))
#    print(i,fcat.GetEntryDescription(i))
#    print(funcGroupList)
#    if len(funcGroupList) > 0:
#        for j in funcGroupList:
#            print(fparams.GetFuncGroup(j).GetProp('_Name'))

#print("Info about functional groups:")
#for i in range(fparams.GetNumFuncGroups()):
#    print(fparams.GetFuncGroup(i).GetProp('_Name'))
    

fgs = identify_functional_groups(m)
print(fgs)