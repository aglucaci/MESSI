#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Given as input

“MACSE2_Out_Codons_recombinants_0.fas.treefile_trimmed_list_RTL_0.8” \
    - a running list of which species is removed after each iteration.

Original alignment = MACSE2_Out_Codons_recombinants_0_renamed.fas



"""


# In[68]:


import sys
import os
from Bio import SeqIO


# In[69]:


#INPUT_ALIGNMENT = "MACSE2_Out_Codons_recombinants_0_renamed.fas"
#RTL_FILE = "MACSE2_Out_Codons_recombinants_0.fas.treefile_trimmed_list_RTL_0.8"

basedir = "p53/"

INPUT_ALIGNMENT = "TP53_refseq_cds_CodonAligned.fas"
RTL_FILE = "PARSED_TREEMER_TRIMMED_LIST.txt"

#RTL_FILE = "TP53_refseq_cds_CodonAligned.nwk_trimmed_list_RTL_0.8" # this is the ones (species) that are left.. not the ones I want.
# I want which species is eliminated (pruned) at each iteration.

#output
path = basedir + "subsampled_fasta"

try:  
    #os.mkdir(path)  
    os.system("mkdir -p " + path)
except OSError as error:  
    print(error)   


# In[70]:


## Helper functioneros -------------------------------------------------------------------------------------------
def subsample_fasta(input_fasta, species_to_exclude, output_fasta):
    
    with open(input_fasta, "r") as handle:
        for n, record in enumerate(SeqIO.parse(handle, "fasta")):
            gene_id = record.id 
            gene_desc = record.description
            gene_seq = str(record.seq).replace("-", "")
            
            if gene_id == species_to_exclude: 
                
                print("# match found")
                continue #skip this one
            else:
                #continue with normal output
                with open(output_fasta, "a") as fh:
                    fh.write(">" + gene_id + "\n" + gene_seq + "\n")
                fh.close()

                
            #end if
        #end for
    #end with
    
                
#end method


# In[71]:



## Main subroutine -----------------------------------------------------------------------------------------------

#Checks to implement, if outputfile_exists, this is bad, means we already ran this script.. delete it.

with open(basedir + RTL_FILE, "r") as fh:
    for n, line in enumerate(fh):
        species = line.strip()
        print("Processing:", n, species)
        
        OUTPUT_FILE = path + "/" + INPUT_ALIGNMENT + "_subsampled_" + str(n) + ".fasta"
        
        
        if n == 0: 
            # first time
            #I want a fasta, with everything in the original fasta, except this species.
            #OUTPUT_FILE = path + "/" + INPUT_ALIGNMENT + "_subsampled_" + str(n) + ".fasta"
            subsample_fasta(basedir + INPUT_ALIGNMENT, species, OUTPUT_FILE)   
        else:
            INPUT = path + "/" + INPUT_ALIGNMENT + "_subsampled_" + str(n-1) + ".fasta"
            #OUTPUT_FILE = path + "/" + INPUT_ALIGNMENT + "_subsampled_" + str(n) + ".fasta"
            
            subsample_fasta(INPUT, species, OUTPUT_FILE) 
        #end if
    
        
    #end for
#end with


## End of cell ----------------------------------------------------------------------------------------------------


# In[ ]:




