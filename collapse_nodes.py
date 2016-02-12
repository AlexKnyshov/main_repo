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
test = tree.as_phyloxml()
tabulate_names(test)
d = {}
for clade in test.find_clades():
	if len(clade.confidences) > 0:
		d[int(clade.name)] = clade.confidences[0].value
count = 0
for dx in sorted(d, reverse=True):
	if d[dx] < 70:
		test.collapse(str(dx))
		count += 1
		print "clade", dx, "collapsed", d[dx]

print count, "nodes collapsed"

#Phylo.draw_ascii(test)
Phylo.write(test, f+".newick.txt", "newick")
