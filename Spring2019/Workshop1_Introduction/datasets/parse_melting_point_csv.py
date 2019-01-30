#!/usr/bin/env python
import pandas as pd
import re

def determine_if_only_CHO(formula):
    keys = formula.keys()
    #print(keys)
    only_CHO = True
    for element in keys:
        if element not in ['C','H','O']:
            only_CHO = False
    return only_CHO

def calc_dou(formula_dict):
    dou = (2.0*float(formula_dict['C'])-float(formula['H']))/2.0
    return dou

#def decode_formula_2(formula_string):
    # Need to be handle strings like
    # C<sub>19</sub>H<sub>15</sub>NO<sub>6</sub>
    

def decode_formula(formula_string):
    #The formula string is the Handbook of Chemistry and Physics style
    # Example C<sub>16</sub>H<sub>20</sub>O<sub>6</sub>P<sub>2</sub>S<sub>3</sub>
    simplified = re.split('</sub>',formula_string)
    #print(simplified)
    split_element = []
    for element in simplified:
        split_element.append(element.split('<sub>'))
        #print(split_element)
    #print(split_element)
    for element in split_element:
        if element[0] == '':
            continue
        else:
            split_uppercase = re.findall('[A-Z][^A-Z]*', element[0])
            if len(split_uppercase) == 1:
                continue
            else:
                for part in split_uppercase[0:len(split_uppercase)-1]: 
                    split_element.append([part,1])
                element[0] = split_uppercase[-1]
    #print(split_element)
    formula = dict()
    for element in split_element:
        #print(element)
        if element[0] == '':
            continue
        else:
            if len(element) > 1:
                if element[0] not in formula.keys():
                    formula[element[0]] = int(element[1])
                else:
                    formula[element[0]] = formula[element[0]] + int(element[1])
            if len(element) == 1:
                if element[0] not in formula.keys():
                    formula[element[0]] = 1
                else:
                    formula[element[0]] = formula[element[0]] + 1
    #print(formula)
    return formula

df = pd.read_csv('melting_points_combined.csv', header=0)

# The headers are mangled: Row is really the name, Synonym is really the formula
# physical form is the melting point, and Mol. form. is the CAS number
#print(df.head(25))
print(df.columns)
#print(df[['Row','Synonym','Mol. form.','Physical form' ]].head(25))
#print(df['Synonym'][1])

col_names=['name','cas','c','h','o','dou','mp']
df2 = pd.DataFrame(columns=col_names)
#print(df2.head())
for i in range(1,len(df['Row'])):
#for i in range(1,5):
    #print(i)
    mf = df['Synonym'][i]
    formula = decode_formula(mf)
    only_CHO = determine_if_only_CHO(formula)
    if only_CHO:
        name = df['Row'][i]
        mp = df['Physical form'][i]
        mp = str(mp).split('(')[0]
        cas = df['Mol. form.'][i]
        c = formula['C']
        h = formula['H']
        if 'O' in formula.keys():
            o = formula['O']
        else:
            o = 0
        dou = calc_dou(formula)
        if mp != 'nan':
            df2 = df2.append({'name':name,
                              'cas':cas,
                              'c':c,
                              'h':h,
                              'o':o,
                              'dou':dou,
                              'mp':mp},ignore_index=True)
        #print(name,cas,c,h,o,dou,mp)
        #tempdf = pd.DataFrame()
        #df2.append(tempdf)
    #print(mf)
    #print(formula)
    #print(only_CHO)
    
    #_ = input()


print(df2)

