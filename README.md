# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).

The TaxonClassifier is used here (folder `RNAs`) to produce the benchmark associated to the paper, but can be applied to any different set of molecules for which some taxonomy information is available. The comparison methods that are considered here are those described in the paper, but any other available method can be included provided that the relative code is available and working.

## Contents

This project contains two Python scripts to process CSV files similar to the ones in folders `RNAs/Distances` and `RNAs/Molecules`. They use agglomerative clustering to reconstruct the taxa associated to a set of molecules according to a distance matrix or pairs of features. The reconstructed taxa are then compared with the ones that are taken from a curated biological taxonomy. 

Both scripts accept as first input a CSV file containing a list o molecules with associated taxonomy information (the taxa of a certain taxonomy rank in a chosen taxonomy). Script `ClusterMatrix.py` accepts as second input a CSV file containing a distance matrix for the molecules in the selected set while the script `ClusterFeatures.py` accept as second input a CSV file containing pairs of features for each molecule in the selected set. In the experiment presented in the folder `RNAs`, the comparison method RAG-2D produces pairs of features while the other methods produce a distance matrix. 

Folder `RNAs/Molecules` contains lists of molecules of the sets Archaea 16S, Archaea 23S, Archaea 5S, Bacteria 16S, Bacteria 23S, Bacteria 5S, Eukaryota 16S, Eukaryota 23S and Eukaryota 5S with associated Phylum Taxa according to the European Nucleotide Archive (ENA) taxonomy <https://www.ebi.ac.uk/ena/browser/home>. The molecule files are available in different formats at <https://doi.org/10.6084/m9.figshare.20731783.v1>. They are part of the Benchmark defined in the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).

Folder `RNAs/Distances` contains the distance matrices or a pair of features computed for the sets Archaea 16S, Archaea 23S, Archaea 5S, Bacteria 16S, Bacteria 23S, Bacteria 5S, Eukaryota 16S, Eukaryota 23S and Eukaryota 5S with the tools Genus, PSMAlign, ASPRAlign, PskOrder and RAG-2D.

Folder `RNAs/Computation` contains scripts and software to run the following comparison methods for pseudoknotted RNA secondary structures: Genus, PSMAlign, ASPRAlign, PskOrder and RAG-2D. They can be run on all the sets Archaea 16S, Archaea 23S, Archaea 5S, Bacteria 16S, Bacteria 23S, Bacteria 5S, Eukaryota 16S, Eukaryota 23S and Eukaryota 5S.

Folder `RNAs/Benchmark` contains the description of the benchmark introduced in the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022). The results of the application of the TaxonClassifier framework to the RNA familes and tools that are considered in the paper are reported and the process to reproduce them is described. The same results are reported and discussed in the paper. 

## Requirements for running the scripts

* Python 3
* Scikit-learn

## Requirements for computing distances/features

* Python 3
* Java java runtime environment (jre) version 8 or higher 
* Bash or a similar shell
* PskOrder binaries (see `Computation/PskOrder/README.md`)
* Matlab, for computing the RAD-2D Fiedler vector features (see `Computation/RAG-2D/README.md`)
* Perl 5.14, for executing PSMAlign (see `Computation/PSMAlign/README.md`)

## Use of the scripts

For example, to obtain the evaluation of the performance of the Genus method, which produces a distance matrix, on 16S Archaea run:

`python3 ClusterMatrix.py RNAs/Molecules/Archaea/16S-Archaea-Phylum.csv RNAs/Distances/Archaea/Genus/16S-genus.csv`

To obtain the evaluation of the performance of the RAG-2D method, which produces pairs of features, on 23S Bacteria run:

`python3 ClusterFeatures.py RNAs/Molecules/Bacteria/23S-Bacteria-Phylum.csv RNAs/Distances/Bacteria/Rag-2d/23S-eigenvalues.csv`

More information about the executions is available at `Benchmark/README.md`

## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.





