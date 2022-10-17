# TaxonClassifier
Evaluation Framework of Comparison Methods that uses Automatic Taxonomy Reconostruction via Agglomerative Clustering

The code and the data available in this project are a support for the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022).

## Benchmark description

We report here the results of the execution of the TaxonClassifier framework on the families of 5S, 16S, 23S rRNAs of Archaea, Bacteria and Eukaryota in the repository at <https://doi.org/10.6084/m9.figshare.20731783.v1>

This benchmark tries to reconstruct the taxa at the **phylum** taxonomy rank according to the European Nucleotide Archive (ENA) taxonomy [1]. The taxonomy was extracted from the SILVA database <https://www.arb-silva.de/>. The phylum taxa of each molecule are reported in the `RNAs/Molecules` subfolder of the main folder. 

These results are also reported in the paper "Automatic Generation of Pseudoknotted RNAs Taxonomy" by Michela Quadrini, Luca Tesei, and Emanuela Merelli (2022) together with a discussion.

## Commands to reproduce the benchmark

Run, from the main directory:

`># python3 ClusterMatrix.py RNAs/Molecules/Kingdom/Family-Phylum.csv RNAs/Distances/Kingdom/Tool/Type-Tool.csv`

For instance:

`># python3 ClusterMatrix.py RNAs/Molecules/Archaea/5S-Archaea-Phylum.csv RNAs/Distances/Archaea/Aspralign/16S-genus.csv`

Produces the following output 

	Method: single
	Rand_score 0.48
	Homogeneity_score 0.16445092468173897
	completeness_score 0.1304306764935216
	Method: complete
	Rand_score 0.48
	Homogeneity_score 0.16445092468173897
	completeness_score 0.1304306764935216
	Method: average
	Rand_score 0.48
	Homogeneity_score 0.16445092468173897
	completeness_score 0.1304306764935216`

For the tool RAG-2D, run, from the main directory:

`># python3 ClusterFeatures.py RNAs/Molecules/Bacteria/23S-Bacteria-Phylum.csv RNAs/Distances/Bacteria/Rag-2d/23S-eigenvalues.csv`

The following output is produced:

	Method: single
	Rand_score 0.583610188261351
	Homogeneity_score 0.2161933823099246
	completeness_score 0.22963103460946435
	Method: complete
	Rand_score 0.5913621262458472
	Homogeneity_score 0.22339298193595716
	completeness_score 0.22007585490161516
	Method: average
	Rand_score 0.5913621262458472
	Homogeneity_score 0.2233929819359571
	completeness_score 0.22007585490161516`


Repeat for all Families: 
- 5S-Archaea
- 16S-Archaea
- 23S-Archaea
- 5S-Bacteria
- 16S-Bacteria
- 23S-Bacteria
- 5S-Eukaryota
- 16S-Eukaryota
- 23S-Eukaryota 

and all the tools: 

- Aspralign
- Genus, 
- Pskorder
- Psmalign
- RAG-2D (eigenvalues)

The files:
- `Benchmarck_Archaea.docx` reports the results on all the Archaea families for all the tools
- `Benchmark_Bacteria.docx` reports the results on all the Bacteria families for all the tools
- `Benchmark_Eukaryota.docx` reports the results on all the Eukaryota families for all the tools
- `Benchmark_All.docx` reports the results on all the families for all the tools

## References

1. Amid, C., Alako, B.T., Balavenkataraman Kadhirvelu, V., Burdett, T., Burgin, J., Fan, J., Harrison, P.W., Holt, S., Hussein, A., Ivanov, E., et al.: The European nucleotide archive in 2019. Nucleic acids research 48(D1), 70–76 (2020)

## Contact Information

Please report any issue to <michela.quadrini@unicam.it> or to Michela Quadrini, Polo Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
