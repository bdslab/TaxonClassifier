#!/bin/bash
# Process a folder containing files generated by DualGraph algorithm of RAG-2D 
# and extract the adjacency matrix in a format that can be easy readable by 
# Matlab. The file is put in the same directory, with the same name and suffix
# _AdjMat.txt
# Usage:
# > ./matrixconverterworkbench.sh <input-folder>

if [[ "$#" -eq 0 ]] ; then
  echo "matrix converter workbench: no input folder given"
  echo "Usage: pseudoorderworkbench.sh <input-folder>"
  exit 1
fi

if ! [[ -d $1 ]] ; then # the folder with the molecules does not exist
  echo "Input folder $1 was not found or is not a folder" 
  exit 1
fi


echo "Processing Molecules in Folder: $1" ...

# Main loop
for mol1 in $(ls $1)
do
    if [[ -d $mol1 ]] ; then
      # do not process directories
      continue
    fi
    
    echo Processing Molecule $mol1 ...
    
    #execute Conversion on the current molecule
    
    java -jar MatrixConverter.jar $1/$mol1 "$1/${mol1}_AdjMat.txt"
done



