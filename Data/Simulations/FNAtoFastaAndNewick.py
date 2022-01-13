#!/usr/bin/env python
# coding: utf-8

# Imports
import sys 
import os
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
from os import path


# In[22]:
DIR = "/home/tuk13147/work/MESSI/Data/Simulations"

onlyfiles = []
onlyfiles = [DIR+"/"+f for f in listdir(DIR) if isfile(join(DIR, f)) and ".replicate" in f and ".newick" not in f]
print("# Found %s files" % str(len(onlyfiles)))


# In[23]:
def parse_newick(filename):
    with open(filename) as f:
        for line in f:
            pass
        #end for
        newick = line # last line
    #end with
    output_newick = filename + ".newick"

    if path.exists(output_newick): 
        return
    else:
        # writeout newick
        text_file = open(output_newick, "w")
        n = text_file.write(newick)
        text_file.close()
    #end if
#end method

for file in onlyfiles:
    #print("# Processing %s" % file)
    parse_newick(file)
    
    #cmd = "head -n-1 filename > newfile"
    output_fasta = file + ".fasta"
    cmd = "head -n-1 " + file + " > " + output_fasta
    #if path.exists(output_fasta): continue
    print(cmd)
    os.system(cmd)
#end for

    


# In[ ]:


#cmd = "head -n-1 filename > newfile"


# In[ ]:


# End of file

