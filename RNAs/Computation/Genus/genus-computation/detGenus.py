#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute the genus of the molecules codified as bpseq

@author: Michela Quadrini

Usage: > python3 detPskOrder.py <dpseq>

"""

from permutation import Permutation
import numpy as np
import os
from glob import glob
import os
import re
import math
import csv
import sys
 
from subprocess import PIPE, Popen
    


if len(sys.argv) < 2:
    print("Usage > python3 detPskOrder.py <bpseq>")
    exit(1)

mol = sys.argv[1]
    
mol2 = re.sub("\_nH.bpseq$", "", os.path.basename(mol))
    
with open(mol, newline='') as csvfile:
         spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
         lista =[]
         a = b = L = unpaired = paired =0
         for row in spamreader:
             L =L+1
             a = row[0]
             b = row[2]
             if(b=='0'): 
                 lista.append(int(a))
                 unpaired = unpaired+1
             else:
                 lista.append(int(b))
                 paired = paired+1
        
         vect_mol = np.array(lista)
                          
pi = Permutation(*vect_mol)
    
lista_sigma = []
for i in range(2,L):
    lista_sigma.append(int(i))
     
lista_sigma.append(1)
vect_sigma = np.array(lista_sigma)
    
sigma = Permutation(*vect_sigma)
    
r =sigma*pi
    
c = len(list(r.to_cycles()))

g = int((paired/2 - c+1)/2)
    
print(g)