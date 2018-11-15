import sys,os  
import csv
import numpy as np

file_in = 'LLNL_GasolineSurr_RedI.c'
with open(file_in, 'r') as infile:
    lines = infile.readlines()
    
idx_beg = 1000000000
idx_end = -1000000000
for i, line in enumerate(lines):
    raw_line = line.strip()
    if "unclassified reactions" in raw_line:
        print line
        idx_beg = i
    if "return" in raw_line and i > idx_beg:
        print line
        idx_end = i
        break

csv_file = 'test.c' 
with open(csv_file, 'w') as outfile:
    for i in range(idx_beg):
        outfile.write(lines[i])

    for i in range(idx_beg+1,idx_end):
        if "q_r[" in lines[i]: 
            new_line = "         q_r = " + lines[i].split("=")[1]
            outfile.write(new_line)
        else:
            outfile.write(lines[i])

    for i in range(idx_end+1,len(lines)):
        outfile.write(lines[i])
