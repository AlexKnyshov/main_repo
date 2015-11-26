import os
import sys
import glob
files = glob.glob("./original_data/T58_L*.phylip")
loci = {}
progbarc = 0
for f in files:
 	print f
 	raxmlf = "-s "+f
 	testcheck = glob.glob("./*.TEST")
 	for test in testcheck:
 		if test.split("/")[1] in os.listdir("./"):
			os.remove(test)
	os.system("raxmlHPC --silent -f x -p 12345 -m GTRGAMMA -n TEST "+raxmlf)
	filename = "RAxML_distances.TEST"
	infile = open(filename, "r")
	lines = []
	for line in infile:
		lines.append(float(line.split("\t")[1]))
	print sum(lines)/len(lines)*100
 	fname = f.split("/")[2]
	locus = fname.split(".")[0]
	loci[locus] = sum(lines)/len(lines)*100
	os.remove("RAxML_distances.TEST")
	os.remove("RAxML_info.TEST")
	os.remove("RAxML_parsimonyTree.TEST")
 	raxmlr = fname+".reduced"
	if raxmlr in os.listdir("./original_data/"):
		os.remove(f+".reduced")
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	hashes = '#' * int(progbar * 0.2)
	spaces = ' ' * (20 - len(hashes))
	print "-------------------------------------"
	print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
	print "-------------------------------------"
	infile.close()
with open("lociGTR.tab", "w") as outfile:
	for lc, lcv in loci.items():
		print >> outfile, lc, "\t", lcv
print "Done"