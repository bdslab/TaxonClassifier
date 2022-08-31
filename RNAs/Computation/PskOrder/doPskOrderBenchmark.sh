#/bin/sh

# Execute the genus pairing for all the folders of the benchmark

python3 PskOrderPairing.py pskorder-computation/Archaea/pskorder16S.csv pskorder-results/Archaea/pskorder16S.csv

python3 PskOrderPairing.py pskorder-computation/Archaea/pskorder5S.csv pskorder-results/Archaea/pskorder5S.csv

python3 PskOrderPairing.py pskorder-computation/Archaea/pskorder23S.csv pskorder-results/Archaea/pskorder23S.csv

python3 PskOrderPairing.py pskorder-computation/Bacteria/pskorder16S.csv pskorder-results/Bacteria/pskorder16S.csv

python3 PskOrderPairing.py pskorder-computation/Bacteria/pskorder5S.csv pskorder-results/Bacteria/pskorder5S.csv

python3 PskOrderPairing.py pskorder-computation/Bacteria/pskorder23S.csv pskorder-results/Bacteria/pskorder23S.csv

python3 PskOrderPairing.py pskorder-computation/Eukaryota/pskorder16S.csv pskorder-results/Eukaryota/pskorder16S.csv

python3 PskOrderPairing.py pskorder-computation/Eukaryota/pskorder5S.csv pskorder-results/Eukaryota/pskorder5S.csv

python3 PskOrderPairing.py pskorder-computation/Eukaryota/pskorder23S.csv pskorder-results/Eukaryota/pskorder23S.csv

