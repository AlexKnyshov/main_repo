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
	#print clade, clade.name
	if len(clade.confidences) > 0:
		#print clade.name
		support = clade.confidences[0].value
		d[int(clade.name)]=support
		#print "support", support
		# if support < 70:
		# 	test.collapse(clade)
		# 	print "clade", clade.name, "collapsed", support
count = 0
for dx in sorted(d, reverse=True):
#	print dx, d[dx]
	if d[dx] < 70:
#		print "low"
		test.collapse(str(dx))
		count += 1
		print "clade", dx, "collapsed", d[dx]

# for clade in test.find_clades():
#  	print clade
#  	if len(clade.confidences) > 0:
#  		support = clade.confidences[0].value
#  		print "support", support

print count

#Phylo.draw_ascii(test)
Phylo.write(test, f+".newick.txt", "newick")

# outer = 0
# inner = 0
# textf = open(f+".newick.txt", "r")
# for l in textf:
# 	for c in l:
# 		print c
# 		if c == "(":
# 			outer += 1
# 		elif c == ")":
# 			inner += 1
# print outer, inner