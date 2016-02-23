from Bio import Phylo
import sys

#-------------------------------------------------------
def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = '%d_%s' % (idx, clade.name)
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names
#-------------------------------------------------------

f = sys.argv[1]

tree = Phylo.read(f, "newick")

tree.collapse_all(lambda c: c.confidence is not None and c.confidence < 70)



#Phylo.draw_ascii(test)
Phylo.write(tree, f+".newick", "newick")
