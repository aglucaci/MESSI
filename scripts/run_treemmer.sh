#!/bin/bash

#run from Treemmer base directory

#1
#python ../Treemmer_v0.3.py testing/MACSE2_Out_Codons_recombinants_0.fas.treefile --verbose 1
#python Treemmer_v0.3.py testing/p53/TP53_refseq_cds_CodonAligned.nwk --verbose 1


#2
#python ../Treemmer_v0.3.py testing/MACSE2_Out_Codons_recombinants_0.fas.treefile --verbose 1 -RTL 0.8
#python Treemmer_v0.3.py testing/p53/TP53_refseq_cds_CodonAligned.nwk --verbose 1 -RTL 0.8 > testing/p53/TP53_treemmer_.txt

#4
python Treemmer_v0.3.py testing/p53/TP53_refseq_cds_CodonAligned.nwk --verbose 2 -RTL 0.8 > testing/p53/TP53_treemmer_verbose2.txt



# end of file
