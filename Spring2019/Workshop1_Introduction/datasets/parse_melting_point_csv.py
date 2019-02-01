#!/usr/bin/env python
import pandas as pd
import re
from rdkit import Chem

import id_fg

df = pd.read_csv('ONSMeltingPoints.tsv',sep='\t')

print(df)
print(df.columns)

