# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).


## Computation of the Genus of RNA secondary structures to determine the Genus distance

The script `genus-computation/detGenus.py` computes the genus of a secondary structure given in BPSEQ format. This script is our implementation of the algorithm presented in Vernizzi et al. (2016) [1]. 

The script `genus-computation/genusworkbench.sh` exectutes the `detGenus.py` script on all the files at a given folder and writes the results in a CSV file.

The script `genus-computation/doBenchmark.sh` executes `genusworkbench.sh` on all the families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota that are present in the repository <https://doi.org/10.6084/m9.figshare.20731783.v1> and reported here in the subfolder `genus-computation/benchmark-bpseq-nH/`. The computed genera for each family are reported in a CSV file of the `genus-computation/genus-results/` subfolder. 

The script `genusPairing.py` pairs all the molecules in a CSV file in order to compute the Genus distance, which is defined as the absolute value of the difference between the two genera. The script `doGenusBenchmark.sh` executes the genusPairing on all the families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota to determing the Genus distances, reported in the subfolder `genus-results/`.

## References

1. Vernizzi, G., Orland, H., Zee, A.: Improved RNA pseudoknots prediction and classification using a new topological invariant. arXiv preprint arXiv:1605.04825 (2016)
2. Giegerich, R., Voß, B., Rehmsmeier, M.: Abstract shapes of RNA. Nucleic acids research 32(16), 4843–4851
(2004)
3. Bon, M., Vernizzi, G., Orland, H., Zee, A.: Topological classification of RNA structures. Journal of molecular
biology 379(4), 900–911 (2008)
4. Reidys, C.M., Wang, R.R.: Shapes of RNA pseudoknot structures. Journal of Computational Biology 17(11),
1575–1590 (2010)


## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
