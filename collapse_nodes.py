from Bio import Phylo
import sys

f = sys.argv[1]

tree = Phylo.read(f, "newick")
Phylo.draw_ascii(tree)