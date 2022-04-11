""" 

MESSI (2022)

Input: multiple sequence alignment (FASTA) 
       and a phylogenetic tree (Newick)

Output: FASTA and Newick files from total number of species (N) to 3, 
         of subsampled alignment and subsampled phylogenetic tree.

N = Total number of species/taxa in the alignment and tree.

Our current procedure will then run selection analyses (MEME) 
    on all subsampled seqs

Author: Alexander G Lucaci

"""

#----------------------------------------------------------------------------
# Imports
#----------------------------------------------------------------------------
import sys
import os
import statistics
from Bio import SeqIO
import glob
import random
import csv
from Bio import SeqIO
import os
from ete3 import Tree
import pandas as pd

#----------------------------------------------------------------------------
# Declares
#----------------------------------------------------------------------------
configfile: 'config.yml'

with open("cluster.json", "r") as in_c:
  cluster = json.load(in_c)
#end with

# Settings from config file
Alignment_file = config["Alignment"]
Tree_file = config["Newick"]
DataSet_Tag = config["Tag"]

Alignment_input = os.path.join("Data", Alignment_file)
Tree_input = os.path.join("Data", Tree_file)

# Get ID's from the fasta file
IDs = []
with open(Alignment_input, "r") as handle:
    for n, record in enumerate(SeqIO.parse(handle, "fasta")):
        IDs.append(record.description)
    #end for
#end with
#print("# IDS:", IDs)

# Report to user
print("# MESSI -- Molecular Evolution SubSampling Investigator ")
print("# @Author: Alexander G Lucaci")
print("# 2022: Version 0.01")

print("# Loading initial FASTA from:", Alignment_input) 
print("# Loading initial NEWICKf from:", Tree_input) 
print("# Total sequences (tips): %i" % len(IDs))
print()

lower_limit = 3 # run until 3 sequences left 

Entries = list(range(lower_limit, len(IDs))) 

print("# Number of subsamples we will examine:", len(Entries))

WD = config["WD"]

OutputDirectory = os.path.join(WD, config["OutputDirectory"], DataSet_Tag)

NP = cluster["__default__"]["ppn"] 

print("# Results will be saved to following folder:", OutputDirectory)

#----------------------------------------------------------------------------
# Helper functions
#----------------------------------------------------------------------------
def output_nwk(tree, excluded_IDs, all_record_IDs, num_seqs, output_dir):
  x = num_seqs - len(excluded_IDs)
  #print("# Output TREE exclusion list length", len(excluded_IDs))
  #output_file = tree + "." + str(count) + ".subsampled.nwk"
  output_file_print = tree + "." + str(x) + ".subsampled.nwk"   
  # Newick, branch pruning logic.
  try:
      t = Tree(tree)
  except:
      t = Tree(tree, format=1)
  #end try
  # Pruning logic here --------------------------------------------
  KeepThese = []
  #for item in excluded_IDs
  for item in all_record_IDs:
      if item not in excluded_IDs:
          KeepThese.append(item)
      #end if
  #end for
  t.prune(tuple(KeepThese)) 
  outfile_tree_print = os.path.join(output_dir, output_file_print.split("/")[-1])
  t.write(format=1, outfile=outfile_tree_print)
  #print("# Saved to (newick file):", outfile_tree_print)
#end method

def output_dna(fasta, excluded_IDs, num_seqs, output_dir):
  x = num_seqs - len(excluded_IDs)
  #print("# Output DNA exclusion list length", len(excluded_IDs))
  output_file = fasta + "." + str(x) + ".subsampled.fasta"
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
  #print("# Saved to (subsampled fasta file):", os.path.join(output_dir, output_file.split("/")[-1]))
#end method

#----------------------------------------------------------------------------
# Rule ALL
#----------------------------------------------------------------------------

rule all:
     input:
         os.path.join(OutputDirectory, Alignment_file + ".dst"),
         os.path.join(OutputDirectory, Alignment_file + ".dst.sorted"),
         expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
         expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries),
         os.path.join(OutputDirectory, Alignment_file + ".MEME.json"),
         expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.MEME.json"), N=Entries)
#end rule all

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# TN93, on codon alignment, can be modified for distance calcs on protein aln
#----------------------------------------------------------------------------
rule TN93:
    input:
       input = Alignment_input
    output:
       output = os.path.join(OutputDirectory, Alignment_file + ".dst")
    shell:
       "tn93 -t 1 -o {output.output} {input.input}"
#end rule

#----------------------------------------------------------------------------
# Sort TN93 Distances
#----------------------------------------------------------------------------
rule Sort_TN93:
    input:
       input = rules.TN93.output.output
    output:
       output = os.path.join(OutputDirectory, Alignment_file + ".dst.sorted")
    run:
       df = pd.read_csv(input.input) 
       df = df.sort_values(by='Distance', ascending=True)
       df.to_csv(output.output)
#end rule

#----------------------------------------------------------------------------
# Subsample
#----------------------------------------------------------------------------
rule Subsample:
    input:
        FASTA     = Alignment_input,
        DISTANCES = rules.Sort_TN93.output.output,
        TREE      = Tree_input
    output:
        fastas =  expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees =  expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    run:
        sorted_distances = pd.read_csv(input.DISTANCES)
        excluded_IDs, count = [], 1
        num_items = sorted_distances.shape[0]
        for index, row in sorted_distances.iterrows():
            # Logic to exit at the end.
            if count == len(IDs):
                #print("# Exiting last iteration, nothing left.")
                break
            #end if
            ID1, ID2, D = row["ID1"], row["ID2"], row["Distance"]
            if ID1 in excluded_IDs and ID2 in excluded_IDs: # both have been previously excluded
                continue
            #end if
            # --- Make a choice ---
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
            # Randomly pick one of them ---
            if to_exclude == "":
                to_exclude = random.choice([ID1, ID2])
            #end if
            # Add to exclusion list ---
            excluded_IDs.append(to_exclude)
            # Save the fasta, and the newick
            output_dna(input.FASTA, excluded_IDs, num_seqs=len(IDs), output_dir=OutputDirectory)
            output_nwk(input.TREE, excluded_IDs, all_record_IDs=IDs, num_seqs=len(IDs), output_dir=OutputDirectory)
            count += 1
        #end for
# end rule

#----------------------------------------------------------------------------
# Selection analyses
#----------------------------------------------------------------------------

rule MEME_FULL: # Run MEME once on the full dataset
    input:
        alignment = Alignment_input,
        tree = Tree_input
    output:
        results = os.path.join(OutputDirectory, Alignment_file + ".MEME.json")
    shell:
        "mpirun -np {NP} HYPHYMPI MEME --alignment {input.alignment} --tree {input.tree} --output {output.results}"
#end rule meme

rule MEME:
    input:
        alignment = os.path.join(OutputDirectory, Alignment_file + ".{Entries}.subsampled.fasta"),
        tree = os.path.join(OutputDirectory, Tree_file + ".{Entries}.subsampled.nwk")   
    output:
        results = os.path.join(OutputDirectory, Alignment_file + ".{Entries}.subsampled.fasta.MEME.json")
    shell:
        "mpirun -np {NP} HYPHYMPI MEME --alignment {input.alignment} --tree {input.tree} --output {output.results}"
#end rule meme

#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
