#!/bin/bash

#Declares
REFERENCE="p53/Human_p53.fasta"


#Main subroutine
for fasta in p53/subsampled_fasta/*.fasta; do
    echo "# Processing: "$fasta

    NUCFILE=$fasta"_nuc.fas"
    PROTFILE=$fasta"_protein.fas"

    # Pre-MSA
    if [ -s $NUCFILE ]; then
        echo "# Nucleotide file exists"
    else
        hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $fasta --reference $REFERENCE --keep-reference No --remove-stop-codons Yes
    fi
    
    prunedNUCFILE=$fasta"_nuc_pruned.fas"
    prunedPROTFILE=$fasta"_protein_pruned.fas"
    
    # Pruned, checking for errors
    if [ -s $prunedNUCFILE ]; then
        echo "# Pruned file exists, checked for - and X's"
    else
        python Prune_nuc_and_protein.py $NUCFILE $PROTFILE
    fi
    
    # MAFFT
    if [ -s $fasta"_protein.msa" ]; then
        echo "# Protein MSA exists"
    else
        #mafft --auto $PROTFILE > $fasta"_protein.msa"
        mafft --auto $prunedPROTFILE > $fasta"_protein.msa"
    fi
    
    # Codon aware alignment
    CODONOUTPUT=$fasta"_codon.alignment"
    
    if [ -s $CODONOUTPUT ]; then
        echo "Codon alignment exists"
    else
        #hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $NUCFILE --output $CODONOUTPUT
        
        echo ""
        echo hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $prunedNUCFILE --output $CODONOUTPUT
        echo ""
        hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $prunedNUCFILE --output $CODONOUTPUT
    fi
    
    # FastTree
    treeOUTPUT=$fasta".nwk"
    if [ -s $treeOUTPUT ]; then
        echo "FastTree exists"
    else
        echo FastTree -gtr -nt $CODONOUTPUT > $treeOUTPUT
        FastTree -gtr -nt $CODONOUTPUT > $treeOUTPUT
    fi
    
    # SA - FEL
    felOUTPUT=$CODONOUTPUT".FEL.json"
    if [ -f $felOUTPUT ]; then
        echo "FEL output exists"
    else
        # Selection Analyses
        #hyphy fel --alignment $CODONOUTPUT --tree $treeOUTPUT --output $felOUTPUT
        echo ""
    fi
    
    # SA - FEL Internal branches
    felOUTPUT=$CODONOUTPUT".FEL_internal.json"
    if [ -f $felOUTPUT ]; then
        echo "FEL internal output exists"
    else
        # Selection Analyses
        #hyphy fel --alignment $CODONOUTPUT --tree $treeOUTPUT --branches Internal --output $felOUTPUT
        echo ""
    fi
    
    # SA - MEME
    MEMEOUTPUT=$CODONOUTPUT".MEME.json"
    if [ -f $MEMEOUTPUT ]; then
        echo "MEME output exists"
    else
        # Selection Analyses
        hyphy MEME --alignment $CODONOUTPUT --tree $treeOUTPUT --output $MEMEOUTPUT
    fi
    
    # SA - MEME Internal branches
    MEMEOUTPUT=$CODONOUTPUT".MEME_internal.json"
    if [ -f $MEMEOUTPUT ]; then
        echo "MEME internal output exists"
    else
        # Selection Analyses
        hyphy MEME --alignment $CODONOUTPUT --tree $treeOUTPUT --branches Internal --output $MEMEOUTPUT
    fi

    # SA - SLAC
    SLACOUTPUT=$CODONOUTPUT".SLAC.json"
    if [ -f $SLACOUTPUT ]; then
        echo "SLAC output exists"
    else
        # Selection Analyses
        echo hyphy SLAC --alignment $CODONOUTPUT --tree $treeOUTPUT --output $SLACOUTPUT
        #hyphy SLAC --alignment $CODONOUTPUT --tree $treeOUTPUT --output $SLACOUTPUT
    fi

done

#Organize outputs
#mkdir -p NEWICKS
#mv p53/subsampled_fasta/*.nwk NEWICKS

mkdir -p p53/FEL
mv p53/subsampled_fasta/*.FEL.json p53/FEL

mkdir -p p53/FEL_internal
mv p53/subsampled_fasta/*.FEL_internal.json p53/FEL_internal

mkdir -p p53/MEME
mv p53/subsampled_fasta/*.MEME.json p53/MEME

mkdir -p p53/MEME_internal
mv p53/subsampled_fasta/*.MEME_internal.json p53/MEME_internal

mkdir -p p53/SLAC
mv p53/subsampled_fasta/*.SLACjson p53/SLAC


# End of file
