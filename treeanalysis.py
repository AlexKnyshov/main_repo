from Bio import Phylo
tree = Phylo.read("./../matrix_052515_02tm1.tre", "newick")
#tree = Phylo.read("./../testtree.tre", "newick")
if tree.is_bifurcating():
	print "bifurcating"
else:
	print "not bifurcating"
counter = 0
ct = 0
test = tree.as_phyloxml()
list4 = test.find_clades()
totconf = 0
lowconf = 0
for l1 in list4:
	counter +=1
	if l1.is_bifurcating():
		ct +=1
	l2 = str(l1._get_confidence())
	if l2.find("Confidence") == 0:
		totconf +=1
		l3=l2.split(",")
		l4 = l3[1][:-1].split("=")
		if int(l4[1]) <70:
			lowconf +=1
print "below 70 per cent", float(lowconf)/(totconf)
print "clades #", counter, ct
#test