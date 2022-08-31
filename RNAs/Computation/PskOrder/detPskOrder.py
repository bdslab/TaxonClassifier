#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scan the input line to find the highest pseudoknot using the following ordering:
    
    #Psk order  Brackets        
    # 0         ()   
    # 1         []   
    # 2         {}   
    # 3         <>  
    # 4         Aa   
    # 5         Bb   
    # 6         Cc 

@author: Luca Tesei and Michela Quadrini

Usage: > python3 detPskOrder.py <dot-bracket-string>

"""

import sys

if len(sys.argv) < 2:
    print("Usage > python3 detPskOrder.py <dot-bracket-string>")
    exit(1)

myString = sys.argv[1]
  
pskOrder = 0

for element in myString:
    if element=='[' and pskOrder < 1:
        pskOrder = 1
    if element=='{' and pskOrder < 2:
        pskOrder = 2
    if element=='<' and pskOrder < 3:
        pskOrder = 3 
    if element=='A' and pskOrder < 4:
        pskOrder = 4
    if element=='B' and pskOrder < 5:
        pskOrder = 5
    if element=='C' and pskOrder < 6:
        pskOrder = 6

print(pskOrder)

        

    