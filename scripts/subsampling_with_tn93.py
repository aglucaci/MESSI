# -*- coding: utf-8 -*-
"""
Subsampling_with_TN93.py
Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/13sFTPkBMccdD_Qpp8DuWhU4UN2ckVT5d
"""
# **Import libraries**"""
import random
import csv
from Bio import SeqIO
import os

"""# **Helper functions**"""

def sort(distances):
  with open(distances, "r") as f:
    csv_input = csv.DictReader(f)
    data = sorted(csv_input, key=lambda row: (row['Distance']))
  return data

def output(fasta, excluded_IDs):
  print("# Output exclusion list length", len(excluded_IDs))
  #output_file = fasta + "_subsampled_" + str(len(excluded_IDs)) + ".fasta"
  output_file = fasta + "." + str(len(excluded_IDs)) + ".subsampled"
    
  # Create empty output file
  with open(os.path.join(output_dir, output_file.split("/")[-1]), 'w') as fp:
    fp.write('')
  #end with
  
  with open(fasta, "r") as handle:
    for n, record in enumerate(SeqIO.parse(handle, "fasta")):
      id = record.id 
      if id in excluded_IDs: 
        continue
      #end if    
      desc = record.description
      seq = record.seq
      #output fasta
      with open(os.path.join(output_dir, output_file.split("/")[-1]), 'a') as f_out:
        SeqIO.write(record, f_out, "fasta")
      #end with
  #end with
  print("# Saved to:", os.path.join(output_dir, output_file.split("/")[-1]))
#end method

def get_num_seqs(fasta):
  count = 0
  with open(fasta, "r") as handle:
    for n, record in enumerate(SeqIO.parse(handle, "fasta")):
        #id = record.id 
        #desc = record.description
        #pseq = record.seq
        count += 1
    #end for
  return count
#end method

"""# **Subsampling logic**"""

def subsample(fasta, distances, num_seqs):
  print("# Sorting")
  sorted_distances = sort(distances) # returns a list of sorted (low to high) tn93 distances
  excluded_IDs = []
  print("# Entering for loop")
  count = 1 # keep track of how many we have removed
  for n, entry in enumerate(sorted_distances):
    
    #if not count < int((num_seqs) * 0.75): # Stop at 25% of of the original number of sequences
    #  break
    #end if
    
    #if not count < 3:
    #  break
    #end if
    
    ID1, ID2, D = entry["ID1"], entry["ID2"], entry["Distance"]
    if ID1 in excluded_IDs and ID2 in excluded_IDs: # both have been previously excluded
      continue
    #end if
    # --- Make a choice
    to_exclude = ""
    if ID1 in excluded_IDs:
      if ID2 not in excluded_IDs:
        to_exclude = ID2
      #end if
    elif ID2 in excluded_IDs:
      if ID1 not in excluded_IDs:
        to_exclude = ID1
      #end if
    #end if
    # If neither IDs is in the exlusion list _ 
    # Randomly pick one of them
    if to_exclude == "":
      to_exclude = random.choice([ID1, ID2])
    #end if
    excluded_IDs.append(to_exclude)
    #print("# Chosen to exclude:", [to_exclude], D, ",".join([ID1, ID2]))
    print("# Iteration:", count)
    print("# Chosen to exclude:", [to_exclude])
    print("# Full entry (for debugging):", ",".join([ID1, ID2, D]))
    count += 1
    output(fasta, excluded_IDs)
    print()
  #end for
#end method

"""# **Main/Driver code**"""
# Test code
#fasta = "Test.aln"
fasta = snakemake.params.fasta
tree = ""
distances = snakemake.params.distances
output_dir = snakemake.params.output_dir

print("## Working on (fasta) --", fasta)
#distances = "Test.dst"
print("## Working on (TN93 Distances) --", distances)

num_seqs = get_num_seqs(fasta)
print("## Number of sequences:", num_seqs)

print("## Starting subsampling")
print()

subsample(fasta, distances, num_seqs)

# End of file

