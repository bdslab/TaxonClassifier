#!/bin/bash
# Process molecules in a folder computing the pseudoknot order and write the results in a csv file
# The default output file, if not given in the command line, is pseudoknotorderworkbenchOutput.csv

# Usage:
# > ./pseudoorderworkbench.sh <input-folder> <output-folder> [<output-file>]

# Put in the following variable the path to the MILP executable file in your machine
# The binaries or the source code of MILP can be obtained from the RNAPDBee group
# http://rnapdbee.cs.put.poznan.pl/ 

pathToMILP="/home/luca/Dropbox/Documenti/Tools/MILP/binaries/milp"

if [[ "$#" -eq 0 ]] ; then
  echo "pseudoorder workbench: no input folder(s) given"
  echo "Usage: pseudoorderworkbench.sh <input-folder> <output-folder> [<output-file>]"
  exit 1
fi

if ! [[ -d $1 ]] ; then # the folder with the molecules does not exist
  echo "Input folder $1 was not found or is not a folder" 
  exit 1
fi

if ! [[ -d $2 ]] ; then # the output folder does not exist
  echo "Output folder $2 was not found or is not a folder" 
  exit 1
fi

# Set output file
outputFile="$2/pseudoknotorderworkbenchOutput.csv"

if [[ -n $3 ]] ; then #the optional outputfile is given 
  outputFile=$3
fi

# Write header line in the csv file

echo "Molecule, Pseudoknot Order, Execution Time [ns]" > $outputFile 

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
    
    mol1Name=${mol1%.bpseq}
        
    # Register starting time
    start=$(date +%s%N) # This does not work on Mac OS
    #start=$(python -c 'import time; print(int(time.time() * 1000000000))')
    
    #execute MILP on the current molecule
    
    $pathToMILP $1/$mol1 > "$2/${mol1Name}_Pseudoorder_Output.txt"
    
    # Register stopping time
    end=$(date +%s%N) # This does not work on Mac OS
    #end=$(python -c 'import time; print(int(time.time() * 1000000000))')
    
    # Compute Execution Time
    executionTime=$(($end-$start))
    
    # Determine the pseudoknot order from the output file
    
    bdExtendedLine=`tail -1 "$2/${mol1Name}_Pseudoorder_Output.txt"`
    
    # Scan the line to find the highest pseudoknot using the following ordering:
    #Psk order  Brackets        
    # 0         ()   
    # 1         []   
    # 2         {}   
    # 3         <>  
    # 4         Aa   
    # 5         Bb   
    # 6         Cc  
    
    pskOrder=`python3 detPskOrder.py $bdExtendedLine`
    
    # Write a line in the csv file
    # Compute the names of the molecules from the name of the file
    
    echo -n "${mol1Name}, " >> $outputFile
    echo -n "${pskOrder}, " >> $outputFile
    echo "$executionTime" >> $outputFile
done



