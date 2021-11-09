""" Description -- Subsampling 2021
Input: alignment and phylogenetic tree
Output: from N to 0, of subsampling alignment, and subsampled phylogenetic tree
N = number of species/taxa in the alignment and tree.

Will then run selection analyses on all subsampled seqs

Currently these include: FitMG94, SLAC, FEL, MEME, FUBAR
Also do BUSTEDS, aBSREL, FMM.. <- (not implemented yet).

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

#----------------------------------------------------------------------------
# Declares
#----------------------------------------------------------------------------
configfile: 'config.yml'

# Settings from config file
Alignment_file = config["Alignment"]
Tree_file = config["Newick"]

Alignment_input = os.path.join("Data", Alignment_file)
Tree_input = os.path.join("Data", Tree_file)

records = list(SeqIO.parse(Alignment_input, "fasta"))

# Report to user
print("# MESSI -- Molecular Evolution SubSampling Investigator ")
print("# @Author: Alexander G Lucaci")
print("# 2021: Version 0.001")

print("# Loading initial FASTA from:", Alignment_input) 
print("# Loading initial NEWICKf from:", Tree_input) 
print("# Total sequences (tips): %i" % len(records))
print()

lower_limit = 3 # run until 3 sequences left ideally, cant go lower.

Entries = list(range(lower_limit, len(records))) 

print("# Number of subsamples we will examine:", len(Entries))

WD = config["WD"]
OutputDirectory = os.path.join(WD, config["OutputDirectory"])


NP = 16

print("# Results will be saved to following folder:", OutputDirectory)

#----------------------------------------------------------------------------
# Helper functions
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Rule ALL
#----------------------------------------------------------------------------
# Do distance calculations

rule all:
     input:
         os.path.join(OutputDirectory, Alignment_file + ".dst"), 
         expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
         expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries),
         expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.FEL.json"), N=Entries)
         #expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.MEME.json"), N=Entries),
         #expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.SLAC.json"), N=Entries)
#end rule all

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# TN93, on codon alignment, can be modified for distance calcs on protein aln
#----------------------------------------------------------------------------
rule tn93:
    input:
       input = Alignment_input
    output:
       output = os.path.join(OutputDirectory, Alignment_file + ".dst")
    conda: 'environment.yml'
    shell:
       "tn93 -t 1 -o {output.output} {input.input}"
#end rule tn93

#----------------------------------------------------------------------------
# Subsampling logic
#----------------------------------------------------------------------------
rule subsampling:
    output:
        fastas =  expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees =  expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    conda: 'environment.yml'
    params:
        fasta = Alignment_input,
        distances = rules.tn93.output.output,
        output_dir = OutputDirectory,
        tree = Tree_input
    script:
        "Scripts/subsampling_with_tn93_fastaANDnewick.py"
# end rule subsampling

#----------------------------------------------------------------------------
# Selection analyses
#----------------------------------------------------------------------------

rule fel:
    input:
        alignment = os.path.join(OutputDirectory, Alignment_file + ".{Entries}.subsampled.fasta"),
        tree = os.path.join(OutputDirectory, Tree_file + ".{Entries}.subsampled.nwk")   
    conda: "environment.yml"
    output:
        results = os.path.join(OutputDirectory, Alignment_file + ".{Entries}.subsampled.fasta.FEL.json")
    shell:
        "mpirun -np {NP} hyphy FEL --alignment {input.alignment} --tree {input.tree} --output {output.results} --ci Yes"
#end rule FEL



"""
rule fel:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.FEL.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            #cmd = "hyphy FEL --alignment" + input.alignments[item] + " --tree " + input.trees[item]
            #cmd = " ".join(["hyphy", "FEL", "--alignment", input.alignments[item], "--tree", input.trees[item], "--ci", "Yes"])
            cmd = " ".join(["mpirun", "-np", "16", "hyphy", "FEL", "--alignment", input.alignments[item], "--tree", input.trees[item], "--ci", "Yes"])
            print(cmd)
            os.system(cmd)        
#end rule FEL
"""
#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
