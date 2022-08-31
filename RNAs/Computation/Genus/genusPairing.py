#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute the genus of the molecules codified as bpseq

@author: Michela Quadrini

Usage: > python3 genusPairing.py <genus-csv-file> <output-csv-file>

"""
import sys
import pandas as pd

if len(sys.argv) < 3:
    print("Usage > python3 genusPairing.py <genus-csv-file> <output-csv-file>")
    exit(1)

# Read the list of molecules
print("Reading", sys.argv[1], "...")
data = pd.read_csv(sys.argv[1], sep=",")
# Initializes the output file
file_distance = open(sys.argv[2], 'a+')
#print(data)
data.columns = ["molecule", "value"]
#print(data.columns)

num_lines = len(data.axes[0])

i=0
while (i<num_lines):
    mol = data["molecule"][i]
    v1 = data["value"][i]
    j=i+1
    while (j<num_lines):
        mol2 = data["molecule"][j]
        v2 = data["value"][j]
        #print("i - j", mol, mol2)
        string_to_add = str(data["molecule"][i]) + str(", ") + str(data["molecule"][j]) + str(", ") +  str(abs(v1-v2))
       
        file_distance.write(string_to_add)
        file_distance.write('\n')
        j=j+1
    
    i=i+1
print(num_lines)
file_distance.close()


    
