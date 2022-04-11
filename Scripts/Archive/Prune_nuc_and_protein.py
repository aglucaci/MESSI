#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Imports

import sys
import os
from Bio import SeqIO


# In[4]:



#NUC_FILE = "p53/TP53_refseq_cds_nuc.fas"
#PROT_FILE = "p53/TP53_refseq_cds_protein.fas"
NUC_FILE = sys.argv[1]
PROT_FILE = sys.argv[2]

print("# pruning:", NUC_FILE)

# In[16]:


nuc_data = []

with open(NUC_FILE, "r") as handle:
        for n, record in enumerate(SeqIO.parse(handle, "fasta")):
            gene_id = record.id 
            gene_desc = record.description
            #gene_seq = str(record.seq).replace("-", "")
            gene_seq = str(record.seq)
            
            if "-" in gene_seq: 
                
                #print(gene_id)
                nuc_data.append(gene_id)
                #Remove the sequences from the nuc.fas, and the protein.fas
                continue
            else:
                #OUTPUT = NUC_FILE.split("/")[-1].replace("_nuc.fas", "_nuc_pruned.fas")
                OUTPUT = NUC_FILE.replace("_nuc.fas", "_nuc_pruned.fas")
             
                with open(OUTPUT, "a") as fh:
                    fh.write(">" + gene_id + "\n" + gene_seq + "\n")
                fh.close()
            #end if
        #end for
#end with
    
            


# In[17]:


prot_data = []

with open(PROT_FILE, "r") as handle:
        for n, record in enumerate(SeqIO.parse(handle, "fasta")):
            gene_id = record.id 
            gene_desc = record.description
            #gene_seq = str(record.seq).replace("-", "")
            gene_seq = str(record.seq)
            
            if "X" in gene_seq: 
                
                #print(gene_id)
                prot_data.append(gene_id)
                #Remove the sequences from the nuc.fas, and the protein.fas
                continue
                
            else:
                #OUTPUT = PROT_FILE.split("/")[-1].replace("_protein.fas", "_protein_pruned.fas")
                OUTPUT = PROT_FILE.replace("_protein.fas", "_protein_pruned.fas")
                with open(OUTPUT, "a") as fh:
                    fh.write(">" + gene_id + "\n" + gene_seq + "\n")
                fh.close()
                
            #end if
        #end for
#end with


# In[15]:
"""

count = len(nuc_data)

running_count = 0

for n, item in enumerate(nuc_data):
    if item in prot_data:
        print(n, item, "# match")
        running_count += 1
    else:
        print(n, item, "# No match")
        
        
    #end if
#end for

print()
print(running_count/count) # check if they match

"""


