#!/usr/bin/env python
import pandas as pd
import re

def decode_formula(formula_string):
    #The formula string is the Handbook of Chemistry and Physics style
    # Example C<sub>16</sub>H<sub>20</sub>O<sub>6</sub>P<sub>2</sub>S<sub>3</sub>
    simplified = re.split('</sub>',formula_string)
    #print(simplified)
    split_element = []
    for element in simplified:
        split_element.append(element.split('<sub>'))
        #print(split_element)
    formula = dict()
    for element in split_element:
        print(element)
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
    print(formula)

df = pd.read_csv('melting_points_1_50.csv', header=0)
# The headers are mangled: Row is really the name, Synonym is really the formula
# physical form is the melting point, and Mol. form. is the CAS number
#print(df.head(25))
print(df.columns)
print(df[['Row','Synonym','Mol. form.','Physical form' ]].head(25))
#print(df['Synonym'][1])
mf = df['Synonym'][1]

df2 = pd.DataFrame()
df2 = df2.append(

decode_formula(mf)

#for compound i
