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

#----------------------------------------------------------------------------
# Declares
#----------------------------------------------------------------------------
configfile: 'config.yml'
Alignment_file = config["Alignment"]
Tree_file = config["Newick"]

Alignment_input = os.path.join("data", Alignment_file)
Tree_input = os.path.join("data", Tree_file)

records = list(SeqIO.parse(Alignment_input, "fasta"))
print("Total sequences: %i" % len(records))
lower_limit = 3 # 3 sequences left.

#Entries = list(range(len(records) - 1, lower_limit))
Entries = list(range(lower_limit, len(records)))
#print(Entries)
OutputDirectory = "results"

print(Entries)
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
         expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.FEL.json"), N=Entries)
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
        "scripts/subsampling_with_tn93_fastaANDnewick.py"
# end rule subsampling

#----------------------------------------------------------------------------
# Selection analyses
#----------------------------------------------------------------------------
rule fel:
    input:
        #alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        #trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
        alignment = rules.subsampling.output.fastas,
        tree = rules.subsampling.output.trees
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.FEL.json"), N=Entries)
    conda: 'environment.yml'
    shell:
        "hyphy FEL --alignment {input.alignment} --tree {input.tree} --output {output.results}"
#end rule FEL




#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
