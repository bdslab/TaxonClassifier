# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).


## Computation of the distances/features for the benchmark

Each folder contains the scripts and the instructions for computing the distances among the molecules or the features  of each molecule

- `ASPRAling`, computation of distances between any pair of molecules in the considered RNA families using the ASPRAlign tool
- `Genus`, computation of the genus of each molecule and the corresponding distance between any pair of molecules in the considered RNA families
- `PskOrder`, computation of the Pseudoknot Order of each molecule and the corresponding distance between any pair of molecules in the considered RNA families
- `PSMAlign`, computation of distances between any pair of molecules in the considered RNA families using the PSMAlign tool
- `RAG-2D`, computation of dual graphs of each molecule in the considered RNA families and computation of the relative Fiedler vector features s and e 


## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
