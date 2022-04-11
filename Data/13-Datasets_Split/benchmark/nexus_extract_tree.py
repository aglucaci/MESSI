import os
import sys
from Bio import Phylo


tree_file = sys.argv[1]
trees = Phylo.parse(sys.argv[1], 'nexus')
for tree in trees:
    #print(tree)
    print("# Saving:", tree.name)
    Phylo.write(tree, tree_file+".newick", "newick")



# End of file
