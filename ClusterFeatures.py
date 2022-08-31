#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TaxonClassifier: execute agglomerative clustering on a set of molecules based on
a pair of features and evaluate the quality of the clustering w.r.t. a known 
clusterization.

Usage: python3 ClusterFeatures.py <molecule-list-csv-file> <eigenvalues-csv-file>

The format of <molecule-list-csv-file> must be <"Id", "Organism", "Taxon"> where
Id is a unique identifier in the file, Organism is a textual description of the 
Id and Taxon is the label associated to Id by a known classification

The format of <eigenvalues-csv-file> must be <"Id", "valueS", "valueS"> where Id
is one Id of the first file and valueS/valueE are features computed for each Id.
The features are used to 

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
    print("Usage: python3 ClusterFeatures.py <molecule-list-csv-file> <eigenvalues-csv-file>")
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

# Create dictionary Id -> Taxon 
taxon_of = dict()
for i in range(len(molecules)) :
    taxon_of[molecules.loc[i].loc['Id']] = molecules.loc[i].loc['Taxon'].strip()

# Create dictionary Id -> Label 
label_of = dict()
for i in range(len(molecules)) :
    label_of[molecules.loc[i].loc['Id']] = molecules.loc[i].loc['Taxon'].strip()

# Read the list of features
print("Reading", sys.argv[2], "...")
distances = pd.read_csv(sys.argv[2], sep=";")

distances.columns = ["Id1", "ValueS", "ValueE"]
print(distances.columns)

# Populate List of Features
ListFeatures =[]
for k in range(len(distances)) :
    i = index_of[distances.loc[k].loc['Id1']]
#    j = index_of[distances.loc[k].loc['Id2']]
    value1 = distances.loc[k].loc['ValueS']
    value2 = distances.loc[k].loc['ValueE']
    ListFeatures.append([value1,value2]) 

# Determine the number of clusters as distinct labels in molecules
n_clusters = len(set(label_of.values()))

# Read the true lables assigned to every Id 
labels_true = list(label_of.values())

# Execute clustering with single linkage and determines the predicted labels for each molecule. In this case features are used and the Euclidean distance is chosen as affinity.

model = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage ='single').fit(ListFeatures)
labels_pred = model.fit_predict(ListFeatures)

# Compute the metrics and print the evaluations
print("Method: single")
#print("Labels Pred", list(labels_pred))
#print("Labels True", list(labels_true))
print("Rand_score", metrics.rand_score(labels_true, labels_pred))
print("Homogeneity_score", metrics.homogeneity_score(labels_true, labels_pred))
print("completeness_score", metrics.completeness_score(labels_true, labels_pred))

# Execute clustering with complete linkage and determines the predicted labels for each molecule. In this case features are used and the Euclidean distance is chosen as affinity.

model = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage ='complete').fit(ListFeatures)
labels_pred = model.fit_predict(ListFeatures)

# Compute the metrics and print the evaluations
print("Method: complete")
#print("Labels Pred", list(labels_pred))
#print("Labels True", list(labels_true))
print("Rand_score", metrics.rand_score(labels_true, labels_pred))
print("Homogeneity_score", metrics.homogeneity_score(labels_true, labels_pred))
print("completeness_score", metrics.completeness_score(labels_true, labels_pred))

# Execute clustering with average linkage and determines the predicted labels for each molecule. In this case features are used and the Euclidean distance is chosen as affinity.

model = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage ='average').fit(ListFeatures)
labels_pred = model.fit_predict(ListFeatures)

# Compute the metrics and print the evaluations
print("Method: average")
#print("Labels Pred", list(labels_pred))
#print("Labels True", list(labels_true))
print("Rand_score", metrics.rand_score(labels_true, labels_pred))
print("Homogeneity_score", metrics.homogeneity_score(labels_true, labels_pred))
print("completeness_score", metrics.completeness_score(labels_true, labels_pred))
