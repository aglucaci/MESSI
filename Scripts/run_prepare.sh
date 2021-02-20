#!/bin/bash

 
## Data retrival  -------------------------------------------------------------------------------------------------------------------------
# Downloaded from https://www.ncbi.nlm.nih.gov/gene/7157/ortholog/?scope=117570&term=TP53
# Download both the transcripts and proteins.
# Use codons.py to get the coding sequence.

## Declares  ------------------------------------------------------------------------------------------------------------------------------
refseqTRANSCRIPTS="p53/TP53_refseq_transcript.fasta"
refseqPROTEINS="p53/TP53_refseq_protein.fasta"

codonsUNALIGNED="p53/TP53_refseq_cds"

proteinMSA="p53/TP53_refseq_cds_protein.msa"

nucleotideUNALIGNED="p53/TP53_refseq_cds_nuc.fas"
proteinUNALIGNED="p53/TP53_refseq_cds_protein.fas"

REFERENCE="p53/Human_p53.fasta"
codons_ALIGNED="p53/TP53_refseq_cds_CodonAligned.fas"


## Main analysis --------------------------------------------------------------------------------------------------------------------------

#0
python p53/codons.py $refseqPROTEINS $refseqTRANSCRIPTS $codonsUNALIGNED
# also manually pruned

#1
echo hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $codonsUNALIGNED --reference $REFERENCE --keep-reference No --remove-stop-codons Yes --skip-realignment Yes
hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $codonsUNALIGNED --reference $REFERENCE --keep-reference No --remove-stop-codons Yes


#2, continue to 3a
#Prune nuc file, remove seqs from protein.fas (Done manual for now)

#3
#mafft --auto $proteinUNALIGNED > $proteinMSA

#3a, continue to 4a
mafft --auto p53/TP53_refseq_cds_protein_pruned.fas > p53/TP53_refseq_cds_protein_pruned.msa

#4
#echo hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $proteinMSA --nucleotide-sequences $nucleotideUNALIGNED --output $codons_ALIGNED
#hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $proteinMSA --nucleotide-sequences $nucleotideUNALIGNED --output $codons_ALIGNED

#4a
hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa p53/TP53_refseq_cds_protein_pruned.msa --nucleotide-sequences p53/TP53_refseq_cds_nuc_pruned.fas --output $codons_ALIGNED

#5
FastTree -gtr -nt $codons_ALIGNED > p53/TP53_refseq_cds_CodonAligned.nwk

#6
hyphy fel --alignment p53/TP53_refseq_cds_CodonAligned.fas --tree p53/TP53_refseq_cds_CodonAligned.nwk

# End of file  ----------------------------------------------------------------------------------------------------------------------------


#for evalue in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do
#    hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $codonsUNALIGNED --reference $REFERENCE --keep-reference No --remove-stop-codons Yes --E $evalue
    #hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $codonsUNALIGNED
    
#    mv $nucleotideUNALIGNED $nucleotideUNALIGNED"_"$evalue
#    mv $proteinUNALIGNED $proteinUNALIGNED"_"$evalue
#done
