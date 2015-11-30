import os
import sys
import glob
import operator
#filepath input
if len(sys.argv) > 1:
	files = glob.glob(sys.argv[1]+"/T58_L*.phylip")
else:
	print "FORMAT: python distances.py [folder with phylip]"
	print "EXAMPLE: python distances.py ./phylip"
	sys.exit()
if len(files) == 0:
	print "no phylip files in the directory"
#starting to process files
loci = {}
progbarc = 0
for f in files:
 	print f
 	raxmlf = "-s "+f
 	#check for raxml leftover files
 	testcheck = glob.glob("./*.TEST")
 	for test in testcheck:
 		if test.split("/")[1] in os.listdir("./"):
			os.remove(test)
	#run RAxML
	os.system("raxmlHPC --silent -f x -p 12345 -m GTRGAMMA -n TEST "+raxmlf)
	#process RAxML output
	filename = "RAxML_distances.TEST"
	infile = open(filename, "r")
	lines = []
	for line in infile:
		lines.append(float(line.split("\t")[1]))
	infile.close()
	print sum(lines)/len(lines)*100
 	fname = f.split("/")[2]
	locus = fname.split(".")[0]
	loci[locus] = sum(lines)/len(lines)*100
	#clean up
	os.remove("RAxML_distances.TEST")
	os.remove("RAxML_info.TEST")
	os.remove("RAxML_parsimonyTree.TEST")
 	raxmlr = fname+".reduced"
	if raxmlr in os.listdir("./original_data/"):
		os.remove(f+".reduced")
	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	hashes = '#' * int(progbar * 0.2)
	spaces = ' ' * (20 - len(hashes))
	print "-------------------------------------"
	print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
	print "-------------------------------------"
#output the final table
with open("lociGTR.tab", "w") as outfile:
	for lc, lcv in sorted(loci.items(), key=operator.itemgetter(1)):
		print >> outfile, lc, "\t", lcv
print "Done"