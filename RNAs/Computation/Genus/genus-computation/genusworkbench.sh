#!/bin/bash
# Process molecules in a folder computing the genus and write the results in a csv file
# The default output file, if not given in the command line, is genusworkbenchOutput.csv

# Usage:
# > ./genusworkbench.sh <input-folder> [<output-file>]

if [[ "$#" -eq 0 ]] ; then
  echo "genus workbench: no input folder(s) given"
  echo "Usage: genusworkbench.sh <input-folder> [<output-file>]"
  exit 1
fi

if ! [[ -d $1 ]] ; then # the folder with the molecules does not exist
  echo "Input folder $1 was not found or is not a folder" 
  exit 1
fi

# Set output file
outputFile=$2

# Write header line in the csv file

echo "Molecule, genus" > $outputFile 

echo "Processing Molecules in Folder: $1" ...

#processedMolecules=""

# Main loop
for mol1 in $(ls $1)
do
    if [[ -d $mol1 ]] ; then
      # do not process directories
      continue
    fi
    
    echo Processing Molecule $mol1 ...
    
    mol1Name=${mol1%_nH.bpseq}
    
    # Register starting time
    # start=$(date +%s%N) # This does not work on Mac OS
    # start=$(python -c 'import time; print(int(time.time() * 1000000000))')
    
        
    # Determine the pseudoknot genus from the output file  
 
    genus=`python3 detGenus.py "$1/$mol1"`
    
    # Register stopping time
    # end=$(date +%s%N) # This does not work on Mac OS
    # end=$(python -c 'import time; print(int(time.time() * 1000000000))')

    # Compute Execution Time
    #executionTime="$(($end-$start))"

    # Write a line in the csv file
    # Compute the names of the molecules from the name of the file
    
    echo -n "${mol1Name}, " >> $outputFile
    echo "${genus}" >> $outputFile
    #echo "$executionTime" >> $outputFile
done



