#!/bin/bash


for file in *.nex; do

    echo $file
    if [ ! -s "$file".fasta ]; then
        seqmagick convert $file "$file".fasta
    else 
        echo "FASTA exists"
    fi
  
    if [ ! -s "$file".newick ]; then
        python nexus_extract_tree.py $file
    else
        echo "NEWICK exists"
    fi

done


exit 0
