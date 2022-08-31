#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TaxonClassifier: execute agglomerative clustering on a set of molecules based on
a distance matrix and evaluate the quality of the clustering w.r.t. a known 
clusterization.

Usage: python3 ClusterMatrix.py <molecule-list-csv-file> <distances-csv-file>

The format of <molecule-list-csv-file> must be <"Id", "Organism", "Taxon"> where
Id is a unique identifier in the file, Organism is a textual description of the 
Id and Taxon is the label associated to Id by a known classification

The format of <distances-csv-file> must be <"Id1", "Id2", "Distance"> where Id1 
and Id2 are two different Ids of the first file and Distance is a floating point
value corresponding to the distance from Id1 and Id2 computed by a chosen 
comparison method. 

The output is given textually as the values of "Rand_score", "Homogeneity_score" 
and "Completeness_score" metrics computed for each executed clustering and for 
each linkage parameter of the clustering algorithm: single, complete and 
average.

@author: Michela Quadrini and Luca Tesei
"""
import math
import sys
import pandas as pd
import numpy as np
from sklearn.metrics import *
from sklearn.cluster import *
from sklearn import metrics
import argparse

if len(sys.argv) != 3 :
    print("Usage: python3 ClusterMatrix.py <molecule-list-csv-file> <distances-csv-file>")
    sys.exit(1)
# Read the list of molecules
molecules = pd.read_csv(sys.argv[1], sep=";")

# Create dictionary Id -> Index
index_of = dict()
for i in range(len(molecules)) :
    index_of[molecules.loc[i].loc['Id']] = i

# Create dictionary Id -> Organism 
organism_of = dict()
for i in range(len(molecules)) :
    organism_of[molecules.loc[i].loc['Id']] = molecules.loc[i].loc['Organism'].strip()

# Create dictionary Id -> Label 
label_of = dict()
for i in range(len(molecules)) :
    label_of[molecules.loc[i].loc['Id']] = molecules.loc[i].loc['Taxon'].strip()

# Read the list of distances
print("Reading", sys.argv[2], "...")
distances = pd.read_csv(sys.argv[2], sep=";")

# Create Distance Matrix
s= (len(molecules),len(molecules))
distance_matrix= np.zeros(s)

# Populate Distance Matrix
for k in range(len(distances)) :
    i = index_of[distances.loc[k].loc['Id1']]
    j = index_of[distances.loc[k].loc['Id2']]
    value = distances.loc[k].loc['Distance']
    distance_matrix[i][j] = value
    distance_matrix[j][i] = value

# Determine the number of clusters as distinct labels in molecules
n_clusters = len(set(label_of.values()))

# Read the true lables assigned to every Id 
labels_true = list(label_of.values())

# Execute clustering with single linkage and determines the predicted labels for each molecule

model = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage ='single').fit(distance_matrix)
labels_pred = model.fit_predict(distance_matrix)

# Compute the metrics and print the evaluations
print("Method: single")
#print("Labels Pred", list(labels_pred))
#print("Labels True", list(labels_true))
print("Rand_score", metrics.rand_score(labels_true, labels_pred))
print("Homogeneity_score", metrics.homogeneity_score(labels_true, labels_pred))
print("Completeness_score", metrics.completeness_score(labels_true, labels_pred))
    
# Execute clustering with complete linkage and determines the predicted labels for each molecule

model = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage ='complete').fit(distance_matrix)
labels_pred = model.fit_predict(distance_matrix)

# Compute the metrics and print the evaluations
print("Method: complete")
#print("Labels Pred", list(labels_pred))
#print("Labels True", list(labels_true))
print("Rand_score", metrics.rand_score(labels_true, labels_pred))
print("Homogeneity_score", metrics.homogeneity_score(labels_true, labels_pred))
print("Completeness_score", metrics.completeness_score(labels_true, labels_pred))


# Execute clustering with average linkage and determines the predicted labels for each molecule

model = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage ='average').fit(distance_matrix)
labels_pred = model.fit_predict(distance_matrix)

# Compute the metrics and print the evaluations
print("Method: average")
#print("Labels Pred", list(labels_pred))
#print("Labels True", list(labels_true))
print("Rand_score", metrics.rand_score(labels_true, labels_pred))
print("Homogeneity_score", metrics.homogeneity_score(labels_true, labels_pred))
print("Completeness_score", metrics.completeness_score(labels_true, labels_pred))
