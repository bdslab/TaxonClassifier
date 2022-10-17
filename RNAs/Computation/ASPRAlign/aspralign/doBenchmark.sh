#! /bin/bash

# Execute ASPRAlign on the molecules of the benchmark

# Comment/Uncomment each line to execute the corresponding family of RNAs

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Archaea/23S -o ../benchmark-results/Archaea/23S-structures.csv \
../benchmark-results/Archaea/23S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Archaea/5S -o ../benchmark-results/Archaea/5S-structures.csv \
../benchmark-results/Archaea/5S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Archaea/16S -o ../benchmark-results/Archaea/16S-structures.csv \
../benchmark-results/Archaea/16S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Bacteria/23S -o ../benchmark-results/Bacteria/23S-structures.csv \
../benchmark-results/Bacteria/23S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Bacteria/5S -o ../benchmark-results/Bacteria/5S-structures.csv \
../benchmark-results/Bacteria/5S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Bacteria/16S -o ../benchmark-results/Bacteria/16S-structures.csv \
../benchmark-results/Bacteria/16S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Eukaryota/23S -o ../benchmark-results/Eukaryota/23S-structures.csv \
../benchmark-results/Eukaryota/23S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Eukaryota/5S -o ../benchmark-results/Eukaryota/5S-structures.csv \
../benchmark-results/Eukaryota/5S-aspralign.csv

java -Xms4G -jar ASPRAlignWorkbench.jar -f \
../benchmark/Eukaryota/16S -o ../benchmark-results/Eukaryota/16S-structures.csv \
../benchmark-results/Eukaryota/16S-aspralign.csv

