#/bin/sh

# Execute the genus pairing for all the folders of the benchmark

python3 genusPairing.py genus-computation/genus-results/Archaea/genus16S.csv genus-results/Archaea/genus16S.csv

python3 genusPairing.py genus-computation/genus-results/Archaea/genus5S.csv genus-results/Archaea/genus5S.csv

python3 genusPairing.py genus-computation/genus-results/Archaea/genus23S.csv genus-results/Archaea/genus23S.csv

python3 genusPairing.py genus-computation/genus-results/Bacteria/genus16S.csv genus-results/Bacteria/genus16S.csv

python3 genusPairing.py genus-computation/genus-results/Bacteria/genus5S.csv genus-results/Bacteria/genus5S.csv

python3 genusPairing.py genus-computation/genus-results/Bacteria/genus23S.csv genus-results/Bacteria/genus23S.csv

python3 genusPairing.py genus-computation/genus-results/Eukaryota/genus16S.csv genus-results/Eukaryota/genus16S.csv

python3 genusPairing.py genus-computation/genus-results/Eukaryota/genus5S.csv genus-results/Eukaryota/genus5S.csv

python3 genusPairing.py genus-computation/genus-results/Eukaryota/genus23S.csv genus-results/Eukaryota/genus23S.csv

