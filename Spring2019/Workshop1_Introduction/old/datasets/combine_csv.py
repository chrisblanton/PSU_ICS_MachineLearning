#!/usr/bin/env python
import glob

list_csv = glob.glob('mp_*.csv')
print(list_csv)

output = open('melting_points_combined.csv','w')

lines = open(list_csv[0],'r').readlines()
output.write(lines[0])

for file in list_csv:
    lines = open(file,'r').readlines()
    for line in lines[1:]:
        output.write(line)
    output.write('\n')
output.close()
