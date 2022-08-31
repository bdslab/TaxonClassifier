#!/bin/bash
# Process pairs of molecules with psmalign and write the results in a csv file
# The default output file, if not given in the command line, is psmalignworkbenchOutput.csv

# Usage:
# > ./psmalignworkbench.sh <folder> [<output-file>]

if [[ "$#" -eq 0 ]] ; then
  echo "psmalign workbench: no input folder given"
  echo "Usage: psmalignworkbench.sh <folder-to-process> [<output-file>]"
  exit 1
fi

if ! [[ -d $1 ]] ; then # the folder with the molecules does not exist
  echo "$1 was not found or is not a folder" 
  exit 1
fi

# Set output file
outputFile="psmalignworkbenchOutput.csv"

if [[ -n $2 ]] ; then #the optional outputfile is given 
  outputFile=$2
fi

# Write header line in the csv file

echo "Molecule 1, Molecule 2, Distance, Execution Time [ns]" > $outputFile 

echo "Processing Molecules in Folder: $1" ...

processedMolecules=""

# Main loop
for mol1 in $(ls $1)
do
  for mol2 in $(ls $1)
  do
    # Processing the pair (mol1, mol2)
    if [[ $mol1 == $mol2 ]] ; then
      # do not process a molecule with itself
      continue
    fi
    
    if [[ $processedMolecules == *"$mol2"* ]] ; then
      # do not process a pair twice
      continue
    fi 

    if [[ -d $mol1 || -d $mol2 ]] ; then
      # do not process directories
      continue
    fi
    
    echo Processing Molecules $mol1 and $mol2 ...
    
    # Register starting time
    start=$(date +%s%N) # This does not work on Mac OS
    #start=$(python -c 'import time; print(int(time.time() * 1000000000))')
    
    #execute psmalign on the current pair of molecules
    
    perl -I. align.pl $1/$mol1 $1/$mol2 > psmalignOutput.txt
    
    # Register stopping time
    end=$(date +%s%N) # This does not work on Mac OS
    #end=$(python -c 'import time; print(int(time.time() * 1000000000))')
    
    # Compute Execution Time
    executionTime=$(($end-$start))
    
    # Determine the alignment cost from the output file
    costLine=$(cat psmalignOutput.txt | grep -E cost)
    cost=${costLine#"Alignment cost: "}
    
    # Write a line in the csv file
    # Compute the names of the molecules from the name of the file
    
    mol1Name=${mol1%_nH.db}
    mol2Name=${mol2%_nH.db}
    echo -n "${mol1Name}, " >> $outputFile
    echo -n "${mol2Name}, " >> $outputFile
    echo -n "${cost}, " >> $outputFile
    echo "$executionTime" >> $outputFile
  done
  # Update processed molecules
  processedMolecules="${processedMolecules}${mol1} "
done

# Remove the output file

rm psmalignOutput.txt


