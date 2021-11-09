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
Alignment_file = config["Alignment"]
Tree_file = config["Newick"]
Alignment_input = os.path.join("Data", Alignment_file)
Tree_input = os.path.join("Data", Tree_file)
records = list(SeqIO.parse(Alignment_input, "fasta"))

# Report to user
print("# MESSI -- Molecular Evolution SubSampling Investigator ")
print("# @Author: Alexander G Lucaci")
print("")

print("# Loading initial FASTA:", Alignment_input) 
print("# Loading initial NEWICK:", Tree_input) 
print("# Total sequences (tips): %i" % len(records))
print()

lower_limit = 3 # run until 3 sequences left ideally, cant go lower..

Entries = list(range(lower_limit, len(records)))
print("# Number of subsamples we will examine:", len(Entries))

OutputDirectory = "Results"
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

rule fel_internal:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.FEL-internal.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            cmd = " ".join(["hyphy", "FEL", "--alignment", input.alignments[item], "--branches", "Internal", "--ci", "Yes", "--tree", input.trees[item]])
            print(cmd)
            os.system(cmd)        
#end rule FEL

rule meme:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.MEME.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            cmd = " ".join(["hyphy", "MEME", "--alignment", input.alignments[item], "--tree", input.trees[item]])
            print(cmd)
            os.system(cmd)
        #end for
#end rule meme

rule slac:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.SLAC.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            cmd = " ".join(["hyphy", "SLAC", "--alignment", input.alignments[item], "--tree", input.trees[item]])
            print(cmd)
            os.system(cmd)
        #end for
#end rule slac

#FitMG="/home/aglucaci/hyphy-analyses/FitMG94/FitMG94.bf"
#NP=8
#echo $HYPHYMP LIBPATH=$RES $FitMG --alignment $FASTA --tree $TREE --rooted No --lrt Yes --type global --frequencies CF3x4 
#$HYPHYMP LIBPATH=$RES $FitMG --alignment $FASTA --tree $TREE --rooted No --lrt Yes --type global --frequencies CF3x4 

FITMG94_bf = "/Users/alex/Documents/hyphy-analyses/FitMG94/FitMG94.bf"

rule fitmg94:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.FITTER.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            cmd = " ".join(["hyphy", FITMG94_bf, "--alignment", input.alignments[item], "--tree", input.trees[item], "--rooted", "No", "--lrt", "Yes", "--type", "global", "--frequencies", "CF3x4"])
            print(cmd)
            os.system(cmd)
        #end for
#end rule fitmg94

rule FUBAR:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.FUBAR.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            #cmd = "hyphy FEL --alignment" + input.alignments[item] + " --tree " + input.trees[item]
            cmd = " ".join(["hyphy", "FUBAR", "--alignment", input.alignments[item], "--tree", input.trees[item]])
            print(cmd)
            os.system(cmd)
        
#end rule fubar

rule BUSTEDS:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.BUSTEDS.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            #cmd = "hyphy FEL --alignment" + input.alignments[item] + " --tree " + input.trees[item]
            cmd = " ".join(["hyphy", "BUSTED", "--alignment", input.alignments[item], "--tree", input.trees[item]])
            print(cmd)
            os.system(cmd)
        
#end rule BUSTEDS

rule aBSREL:
    input:
        alignments = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta"), N=Entries),
        trees = expand(os.path.join(OutputDirectory, Tree_file + ".{N}.subsampled.nwk"), N=Entries)
    output:
        results = expand(os.path.join(OutputDirectory, Alignment_file + ".{N}.subsampled.fasta.aBSREL.json"), N=Entries)
    run:
        for item in range(len(input.alignments)):
            print(input.alignments[item], input.trees[item])
            #cmd = "hyphy FEL --alignment" + input.alignments[item] + " --tree " + input.trees[item]
            cmd = " ".join(["hyphy", "ABSREL", "--alignment", input.alignments[item], "--tree", input.trees[item]])
            print(cmd)
            os.system(cmd)
        
#end rule aBSREL

#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
