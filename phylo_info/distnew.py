import os
import sys
import glob
import operator
import socket
from Bio import AlignIO
#filepath input
if len(sys.argv) == 4:
	ftype = sys.argv[1]
	inputfolder = sys.argv[3]
	files = glob.glob(inputfolder+"/*")
	mode = sys.argv[2]
else:
	print "FORMAT: python distances.py [ftype: -fas, -phy] [mode: -p, -115, -gtr] [folder with phylip]"
	print "EXAMPLE: python distances.py -phy -gtr ./phylip"
	sys.exit()
if len(files) == 0:
	print "no phylip files in the directory"

#starting to process files
loci = {}
progbarc = 0

if ftype == "-phy":
	format = "phylip-relaxed"
elif ftype == "-fas":
	format = "fasta"

if mode == "-p" or mode == "-115":
	for f in files:
		print f
		infile = open(f, "r")
		seqs = {}
		for seq in AlignIO.read(infile, format):
		    seqs[seq.id] = str(seq.seq).upper()
		#pairwise
		names =[]
		for key in seqs.keys():
			names.append(key)
		distlist = []
		for num in range(len(names)-1):
			pair1 = []
			pair1 += seqs[names[num]]
			for p2 in range(num+1, len(names)):
				pair2 = []
				pair2 += seqs[names[p2]]
				counter = 0
				countert = 0
				counterindel = 0
				startpos = 0
				endpos = len(pair1)
				#forward trim
				needskip = True
				for index in range(len(pair1)):
					if pair1[index] != "-" and pair2[index] != "-" and pair1[index] != "N" and pair2[index] != "N":
						startpos = index
						needskip = False
						break
				#reverse trim
				for index in range(len(pair1)-1, -1, -1):
					if pair1[index] != "-" and pair2[index] != "-" and pair1[index] != "N" and pair2[index] != "N":
						endpos = index
						needskip = False
						break
				if needskip:
					val = "N/A"
				else:
					#checking the differences
					for index in range(startpos, endpos+1):
						if pair1[index] == "A" and pair2[index] == "G":
							countert += 1
						elif pair1[index] == "G" and pair2[index] == "A":
							countert += 1
						elif pair1[index] == "A" and pair2[index] == "C":
							counter +=1
						elif pair1[index] == "A" and pair2[index] == "T":
							counter +=1
						elif pair1[index] == "T" and pair2[index] == "A":
							counter +=1
						elif pair1[index] == "T" and pair2[index] == "G":
							counter +=1
						elif pair1[index] == "C" and pair2[index] == "T":
							countert +=1
						elif pair1[index] == "T" and pair2[index] == "C":
							countert +=1
						elif pair1[index] == "C" and pair2[index] == "A":
							counter +=1
						elif pair1[index] == "C" and pair2[index] == "G":
							counter +=1
						elif pair1[index] == "G" and pair2[index] == "C":
							counter +=1
						elif pair1[index] == "G" and pair2[index] == "T":
							counter +=1
						elif pair1[index] == "-" and pair2[index] != "-":
							counterindel +=1
						elif pair2[index] == "-" and pair1[index] != "-":
							counterindel +=1
					#calculating the distance
					if mode == "-115":
						val = round((float(counter*5+countert+counterindel)/(endpos-startpos+1)*100), 2) #1-1-5 costmatrix, percent
					else:
						val = (float(counter+countert+counterindel)/(endpos-startpos+1)) #normal p-distance as in R package {ape} with pairwise deletion except weird indel counter
				if val != "N/A":
					distlist.append(val)
		print "average distance", sum(distlist)/len(distlist)
		fname = f.split("/")
		locus = fname[2].split(".")[0]
		loci[locus] = sum(distlist)/len(distlist)
		infile.close()
		#progress bar
		progbarc +=1
		progbar = int(round(float(progbarc)/len(files)*100, 0))
		hashes = '#' * int(progbar * 0.2)
		spaces = ' ' * (20 - len(hashes))
		print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
elif mode == "-gtr":
	for f in files:
	 	print f
	 	raxmlf = "-s "+f
	 	#check for raxml leftover files
	 	testcheck = glob.glob("./*.TEST")
	 	for test in testcheck:
	 		if test.split("/")[1] in os.listdir("./"):
				os.remove(test)
		#run RAxML
		hostname = socket.gethostname()
		if hostname == "pigeon":
			execfile('/usr/share/Modules/init/python.py')
			module('load','RAxML/8.2.3')
			os.system("cd ~/gen220project")
			os.system("raxmlHPC-PTHREADS-SSE3 -T 16 --silent -f x -p 12345 -m GTRGAMMA -n TEST "+raxmlf)
		else:
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
		if raxmlr in os.listdir(inputfolder+"/"):
			os.remove(f+".reduced")
		#progress bar
		progbarc +=1
		progbar = int(round(float(progbarc)/len(files)*100, 0))
		hashes = '#' * int(progbar * 0.2)
		spaces = ' ' * (20 - len(hashes))
		print "-------------------------------------"
		print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
		print "-------------------------------------"
if mode == "-p":
	outf = "lociraw.tab"
elif mode == "-115":
	outf = "loci115.tab"
elif mode == "-gtr":
	outf = "locigtrn.tab"
with open(outf, "w") as outfile:
	for lc, lcv in sorted(loci.items(), key=operator.itemgetter(1)):
		print >> outfile, lc, "\t", lcv
print "Done"