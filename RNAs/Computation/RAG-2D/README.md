# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).


## Use of the RAG-2D (RNA-AS-GRAPHS 2D) approach for computing features

Schlick's Lab <http://www.biomath.nyu.edu/?q=rag/home> uses dual graphs to formalize RNA secondary structures with pseudoknots [1, 2, 3] within the RNA As Graph approach. Dual graphs describe the connectivity of an RNA secondary structure without specifying the geometric aspects. They use the Laplacian matrix of the connectivity graph of the structure. In [4] Zhu and Schlick used Fiedler vectors, which are the eigenvectors corresponding to the second smallest eigenvalues of the Laplacain matrices, to define two features, namely s and e, that better reflect the graph topology. They can be used to define a distance between two structures. 

To compute the distances, first the dual graphs of all the molecules must be determined. Folder `DualGraphs` contains the scripts to do so. As courtesy of the Schlick's Lab we obtained the script `DualGraphs/dualGraphs.py`, a variant of the original Python script in <https://github.com/Schlicklab/DualGraphCodes> that permits to process all the files in a given folder instead of one input file (which is the functionality of the original code). The script `DualGraphs/doBenchmark.sh` executes `dualGraphs.py` on all the families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota that are present in the repository <https://doi.org/10.6084/m9.figshare.20731783.v1> and reported here in the subfolder `DualGraphs/benchmark/` taking the ct format without header. For each molecule in each family the dual graph is computed and reported in the corresponding file in the subfolder `benchmark-results/`. 

The output contained too much information for the subsequent step of the computation of the distance. We implemented a little converter in Java (`DualGraphs/MatrixConverter.jar`) that takes a file containing the dual graph and converts it to a simplified version that can be later used as input in the computation of the s and e features of each molecule. The script `DualGraphs/matrixconverterworkbench.sh` applies the converter to all the dual graphs in in the subfolder `benchmark-results/`. 

The second step is to determine the features s and e for each molecule using the converted dual graph. Folder `FeaturesComputation` contains the code to do so. The converted dual graphs computed in the previous step are copied in the subfoleder `FeaturesComputation/benchmark-adjmatrices`. The Matlab script `FeaturesComputation/FeaturesComputationDoBenchmark.m`, courtesy of the Schlick's Lab and partially modified by us, computes the features e and s for each molecule and stores them in CSV files in `FeaturesComputation/benchmark-results` subfolder. 

These features can then be used directly to clusterise the molecule using the `
ClusterFeatures.py` script of the TaxonClassifier.

## References

1. Gan, H.H., Fera, D., Zorn, J., Shiffeldrim, N., Tang, M., Laserson, U., Kim, N., Schlick, T.: RAG: RNA-As-Graphs database—concepts, analysis, and features. Bioinformatics 20(8), 1285–1291 (2004)
2. Gan, H.H., Pasquali, S., Schlick, T.: Exploring the repertoire of RNA secondary motifs using graph theory; implications for RNA design. Nucleic Acids Research 31(11), 2926–2943 (2003)
3. Fera, D., Kim, N., Shiffeldrim, N., Zorn, J., Laserson, U., Gan, H.H., Schlick, T.: RAG: RNA-As-Graphs web resource. BMC bioinformatics 5(1), 1–9 (2004)
4. Zhu, Q., Schlick, T.: A Fiedler Vector Scoring Approach for Novel RNA Motif Selection. The Journal of Physical Chemistry B 125(4), 1144–1155 (2021)


## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
