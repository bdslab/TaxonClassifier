# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).


## Use of the tool ASPRAlign for computing distances

The tool [ASPRAlign](https://github.com/bdslab/aspralign) can align RNA Secondary Structures with arbitrary pseudoknots. The cost of the best alignment between two structures is called ASPRA distance.

The script `aspralign/dobenchmark.sh` executes ASPRAlignWorkbench on all the families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota that are present in the repository <https://doi.org/10.6084/m9.figshare.20731783.v1> and reported here in the subfolder `benchmark/` taking the dot-bracket-letter format without header. For each pair of molecules in each family the ASPRA distance is computed and reported in a CSV file of the `benchmark-results/` subfolder. 

## References

1. Quadrini, M., Tesei, L., Merelli, E.: An algebraic language for RNA pseudoknots comparison. BMC Bioinformatics 20(4), 161 (2019).
2. Quadrini, M., Tesei, L., Merelli, E.: ASPRAlign: a tool for the alignment of RNA secondary structures with arbitrary pseudoknots. Bioinformatics 36(11), 3578–3579 (2020).


## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
