#!/bin/bash
clear


## Declares

BASEDIR=/Users/user/Documents/Subsampling
DATA=$BASEDIR/data/TP53_refseq_cds_CodonAligned.nwk
OUTPUT_DIR=$BASEDIR/analysis
mkdir -p $OUTPUT_DIR

## Software
PYTHON=/Users/user/opt/anaconda3/bin/python
TREEMMER=../Treemmer/Treemmer_v0.3.py


## Main
# Creates LD and Decay plot pdf.
#echo $PYTHON $TREEMMER $DATA
#$PYTHON $TREEMMER $DATA

# The output of this is the trimmed list, and trimmed tree, basically what remains at the 80% RTL.
#echo $PYTHON $TREEMMER $DATA -RTL 0.8
#$PYTHON $TREEMMER $DATA -RTL 0.8

# Default here...
#prints tree at each iteration to screen (in ascii form).
echo $PYTHON $TREEMMER $DATA -RTL 0.8 --verbose 2 -c 2 -pc > "$DATA"_verbose2.txt
$PYTHON $TREEMMER $DATA -RTL 0.8 --verbose 2 -c 2 -pc > "$DATA"_verbose2.txt

## Move files
#echo mv $DATA"_res*" $OUTPUT_DIR/
#echo "mv $DATA_trimmed_list_* $OUTPUT_DIR/"
#mv $DATA_trimmed_list_* $OUTPUT_DIR/


#mv $DATA_* $OUTPUT_DIR

for file in "$DATA"_*; do
  echo $file
  mv $file $OUTPUT_DIR
done


# END OF FILE
