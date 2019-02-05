#!/usr/bin/env python
import pandas as pd
import re
from rdkit import Chem

import id_fg

def get_chem_form_inchi(inchi):
    #InChI are strings like this
    #InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3
    # That is ethanol
    #All InChI start with InchI=
    # The next part is the version and standard number
    # After the / is the chemical formula
    # After the next next /c is the connectivity of the heavy atoms
    # Then after /h is the is the hydrogren connectivity.
    # There are additional layers as needed.
    inchi_string =  inchi.split('=')[-1]
    chem_form = inchi_string.split('/')[1]
    return chem_form

def split_chem_form(chem_form):
    set = re.findall(r'([A-Z][a-z]*)(\d*)', chem_form)
    form_dict = {}
    for element in set:
        if element[0] in form_dict.keys():
            count = '1'
            if element[1] != '':
                form_dict[element[0]] = int(form_dict[element[0]]) + int(element[1])
            else:
                form_dict[element[0]] = int(form_dict[element[0]]) + int(count)
        if element[0] not in form_dict.keys():
            count = '1'
            if element[1] != '':
                form_dict[element[0]] = element[1]
            else:
                form_dict[element[0]] = count
    return form_dict
    

df = pd.read_csv('ONSMeltingPoints.tsv',sep='\t')

print(df)
print(df.columns)

print(df[['Ave','SMILES','name']].head())
for i in range(len(df['SMILES'])):
    m = Chem.MolFromSmiles(df['SMILES'][i])
    inchi = Chem.MolToInchi(m)
    chem_form = get_chem_form_inchi(inchi)
    form_dict = split_chem_form(chem_form)
    print(df['SMILES'][i], inchi, chem_form, form_dict)
    #inp = input('Continue [Y/n]')
    #if inp == 'n':
    #    break