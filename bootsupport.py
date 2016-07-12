import sys
import glob
from Bio import Phylo
from Bio.Phylo.Consensus import *
print "reading the best tree..."
tree = Phylo.read(sys.argv[1], "newick")
Phylo.draw_ascii(tree)
print "reading the bootstrap trees..."
files = glob.glob("./*.r?_nhPhymlEq")
treelist = []
for f in files:
	bstree = Phylo.read(f, 'newick')
	treelist.append(bstree)
	Phylo.draw_ascii(bstree)
print "calculating the support..."
support_tree = get_support(tree, treelist)
Phylo.draw_ascii(support_tree)
#print support_tree
print "writing the file..."
Phylo.draw(support_tree)
#Phylo.draw_graphviz(tree)
#pylab.show()
#pylab.savefig('apaf.png')
#Phylo.write(support_tree, "BStree.tre", "newick")
print "done"


