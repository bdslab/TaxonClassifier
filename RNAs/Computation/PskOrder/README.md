# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).


## Use of the Pseudoknot Order (PskOrder) for computing distances

The Pseudoknot Order measures the structural complexity of a secondary structure with pseudoknots. It can be computed with several algorithms. With the courtesy of the Antczak et al. Lab (see <http://rnapdbee.cs.put.poznan.pl/>) we obtained the executable code of the MILP method implementation [3].


The script `pseudoorderworkbench.sh` executes the MILP mehtod on a folder of input molecules and writes the ouptut of the tool in a specified output folder. At the same time it determines the PskOrder of each molecule and writes it on the specified output CSV file. To use the script the user must indicate the path to the MILP executable file in his/her machine in the variable `pathToMILP` inside the script. This can be obtained from the Antczak et al. Lab.

The script `doBenchmark.sh` executes `pseudoorderworkbench.sh` on all the families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota that are present in the repository <https://doi.org/10.6084/m9.figshare.20731783.v1>. The outpu is written on subfolder `pskorder-computation`.

The script `PskOrderPairing.py` pairs all the molecules in a CSV file in order to compute the PskOrder distance, which is defined as the absolute value of the difference between the two orders. The script `doPskOrderBenchmark.sh` executes the `PskOrderPairing` on all the families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota to determing the PskOrder distances, reported in the subfolder `pskorder-results/`.

## References

1. Antczak, M., Zok, T., Popenda, M., Lukasiak, P., Adamiak, R.W., Blazewicz, J., Szachniuk, M.: RNApdbee—a webserver to derive secondary structures from pdb files of knotted and unknotted RNAs. Nucleic acids research 42(W1), 368–372 (2014)
2. Antczak, M., Popenda, M., Zok, T., Zurkowski, M., Adamiak, R.W., Szachniuk, M.: New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34(8), 1304–1312 (2018)
3. Zok, T., Badura, J., Swat, S., Figurski, K., Popenda, M., Antczak, M.: New models and algorithms for RNA pseudoknot order assignment. International Journal of Applied Mathematics and Computer Science 30(2), 315–324 (2020)


## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
