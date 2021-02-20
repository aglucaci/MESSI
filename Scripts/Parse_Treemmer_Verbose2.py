#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys


# In[2]:


TREEMMERFILE="p53/TP53_treemmer_verbose2.txt"


# In[3]:


species_list = []


with open(TREEMMERFILE) as fh:
    species_next = False
    
    for n,line in enumerate(fh):
        #print(n, [line])
        
        if species_next == True:
            #keep this line.
            species_list.append(line.split("\t")[0]) #there some more information I dont need, BL?
            species_next = False
        
        if "leaf to prune:" in line:
            species_next = True
            
            
        


# In[4]:


species_list


# In[5]:


len(species_list)


# In[6]:


with open("p53/PARSED_TREEMER_TRIMMED_LIST.txt", "w") as handle:
    for item in species_list:
        handle.write(item + "\n")
    


# In[ ]:




