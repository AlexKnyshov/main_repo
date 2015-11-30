from Bio import Phylo
import os
import sys
import glob
import getpass
#filepath input
if len(sys.argv) == 3:
	files = glob.glob(sys.argv[1]+"/*.phylip")
	#files = ['./original_data/T58_L1.phylip']
	bootreps = " -N "+sys.argv[2]
else:
	print "FORMAT: python treeanalysis.py [folder with phylip] [bootstrap replicates]"
	print "EXAMPLE: python treeanalysis.py ./phylip 1000"
	sys.exit()
if len(files) == 0:
	print "no phylip files in the directory"
#files = ['./original_data/T58_L1.phylip']
#starting to process files
progbarc = 0
result = {}

for f in files:
	print f
 	raxmlf = " -s "+f
 	#check for raxml leftover files
 	testcheck = glob.glob("./*.TEST")
 	for test in testcheck:
 		if test.split("/")[1] in os.listdir("./"):
			os.remove(test)
	#run RAxML
	user = getpass.getuser()
	if user == "aknys001":
		os.system("module load RAxML/8.2.3")
		os.system("cd ~/gen220project")
		os.system("raxmlHPC-HYBRID --silent -f a -c 25 -p 12345 -x 12345 -m GTRCAT -n TEST"+raxmlf+bootreps)
	else:
		os.system("raxmlHPC --silent -f a -c 25 -p 12345 -x 12345 -m GTRCAT -n TEST"+raxmlf+bootreps)
	tree = Phylo.read("RAxML_bipartitions.TEST", "newick")
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