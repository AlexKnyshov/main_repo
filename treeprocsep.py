from Bio import Phylo
import os
import sys
import glob
if len(sys.argv) == 2:
	files = glob.glob(sys.argv[1]+"/*.tre")
else:
	print "FORMAT: python treeanalysis.py [folder with trees]"
	print "EXAMPLE: python treeanalysis.py ./trees"
	sys.exit()
if len(files) == 0:
	print "no trees in the directory"
#starting to process files
progbarc = 0
result = {}
#test
for f in files:
	tree = Phylo.read(f, "newick")
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
	fname = f.split("/")[2]
	tr = fname.split(".")[0]
	result[tr] = float(lowconf)/(totconf)
	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	hashes = '#' * int(progbar * 0.2)
	spaces = ' ' * (20 - len(hashes))
	print "-------------------------------------"
	print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
	print "-------------------------------------"
#output the final table
with open("treesfilter.tab", "w") as outfile:
	for tc, tcv in result.items():
		print >> outfile, tc, "\t", tcv
print "Done"