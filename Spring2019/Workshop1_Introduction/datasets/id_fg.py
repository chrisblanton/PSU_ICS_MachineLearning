#!/usr/bin/env python
from __future__ import print_function
from rdkit import Chem

# Useful references
# Ruggeri and Takahama has a very useful table of SMARTS
# Technical Note: Development of chemoinformatic tools to enumerate
#functional groups in molecules for organic aerosol characterization
#Giulia Ruggeri and Satoshi Takahama
# https://core.ac.uk/download/pdf/145658301.pdf
#
# DAYLIGHT tutorials
# http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
# http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

# ROH
def is_alcohol(mol):
    # alcohol, not in carboxylic acid
    #alcohol_smarts = '[$([CX4]O)]'
    alcohol_smarts = '[C,c;!$(C=O)][OX2H]'
    alcohol = Chem.MolFromSmarts(alcohol_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,alcohol)) > 0:
        value = True
    return value

# RCOOH
def is_cooh(mol):
    cooh_smarts = '[CX3](=O)[OX2H1]'
    cooh = Chem.MolFromSmarts(cooh_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,cooh)) > 0:
        value = True
    return value

#R(C=O)R
def is_ketone(mol):
    ketone_smarts = '[#6][CX3](=O)[#6]'
    ketone = Chem.MolFromSmarts(ketone_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,ketone)) > 0:
        value = True
    return value

# ROR
def is_ether(mol):
    #ether_smarts = '[OD2]([#6])[#6]'
    ether_smarts = '[OD2]([#6;!$(C=O)])[#6;!$(C=O)]'
    ether = Chem.MolFromSmarts(ether_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,ether)) > 0:
        value = True
    return value

# RCOOR
# This will not hit anhydrides
def is_ester(mol):
    #ester_smarts = '[#6][CX3](=O)[OX2H0][#6]'
    ester_smarts = '[CX3H1,CX3](=O)[OX2H0][#6;!$([C]=[O])]'
    ester = Chem.MolFromSmarts(ester_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,ester)) > 0:
        value = True
    return value

# R(C=O)O(C=O)R
def is_anhydride(mol):
    anh_smarts = '[CX3](=O)[O][CX3](=O)'
    anh = Chem.MolFromSmarts(anh_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,anh)) > 0:
        value = True
    return value

# RCHO, But does not catch formaldehyde
def is_aldehyde(mol):
    ald_smarts = '[CX3H1](=O)[#6]'
    ald = Chem.MolFromSmarts(ald_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,ald)) > 0:
        value = True
    return value

# Aromatic
def is_aromatic(mol):
    ar_smarts = 'c'
    ar = Chem.MolFromSmarts(ar_smarts)
    value = False
    if len(Chem.Mol.GetSubstructMatch(mol,ar)) > 0:
        value = True
    return value

def find_longest_alphatic_chain(mol):
    alp_smarts = '[C]'
    alp = Chem.MolFromSmarts(alp_smarts)
    value = Chem.Mol.GetSubstructMatches(mol,alp)
    return value
    




        