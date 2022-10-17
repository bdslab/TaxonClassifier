# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).


## Use of the tool PSMAlign for computing distances

Progressive Stem Matching Alignment (PSMAlign) [1] is an alignment method that efficiently computes distances between RNA secondary structures with arbitrary pseudoknots. A corresponding tool PSMAlign was developed and is available at <http://homepage.cs.latrobe.edu.au/ypchen/psmalign/>.


To compute the distances, first download the PSMAlign tool from <http://homepage.cs.latrobe.edu.au/ypchen/psmalign/psmalign.tar.gz> and extract it in the folder `psmalign`. 

The script `psmalign/psmalignworkbench.sh` executes PSMAlign on all pairs of molecules that are present in a given input folder. The molecules must be in dot-bracket-letter format without header. The resulting distances are written in a CSV output file. 

The script `psmalign/dobenchmark.sh` executes `psmalignworkbench.sh` on all the
families of 16S, 5S and 23S RNAs of Archaea, Bacteria and Eukaryota that are present in the repository <https://doi.org/10.6084/m9.figshare.20731783.v1> and reported here in the folder `benchmark/`. The folder `benchmark-results/` contains the CSV files with all the computed distances.

## References

1. Chiu, J.K.H., Chen, Y.-P.P.: Pairwise RNA secondary structure alignment with conserved stem pattern. Bioinformatics 31(24), 3914–3921 (2015)


## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
