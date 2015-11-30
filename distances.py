from Bio import SeqIO
import glob
import sys
import os
import operator
#filepath input
if len(sys.argv) > 1:
	files = glob.glob(sys.argv[1]+"/*.fas")
else:
	print "FORMAT: python distances.py [folder with fasta]"
	print "EXAMPLE: python distances.py ./fasta"
	sys.exit()
if len(files) == 0:
	print "no fasta files in the directory"
#starting to process files
loci = {}
progbarc = 0
for f in files:
	print f
	infile = open(f, "r")
	seqs = {}
	for seq in SeqIO.parse(infile, "fasta"):
	    seqs[seq.id] = str(seq.seq)
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
						counter +=5
					elif pair1[index] == "A" and pair2[index] == "T":
						counter +=5
					elif pair1[index] == "T" and pair2[index] == "A":
						counter +=5
					elif pair1[index] == "T" and pair2[index] == "G":
						counter +=5
					elif pair1[index] == "C" and pair2[index] == "T":
						countert +=1
					elif pair1[index] == "T" and pair2[index] == "C":
						countert +=1
					elif pair1[index] == "C" and pair2[index] == "A":
						counter +=5
					elif pair1[index] == "C" and pair2[index] == "G":
						counter +=5
					elif pair1[index] == "G" and pair2[index] == "C":
						counter +=5
					elif pair1[index] == "G" and pair2[index] == "T":
						counter +=5
					elif pair1[index] == "-" and pair2[index] != "-":
						counterindel +=1
					elif pair2[index] == "-" and pair1[index] != "-":
						counterindel +=1
				#calculating the distance
				val = round((float(counter+countert+counterindel)/(endpos-startpos+1)*100), 2) #1-1-5 costmatrix, percent
				#val = (float(counter+countert+counterindel)/(endpos-startpos+1)) #normal p-distance as in R package {ape} with pairwise deletion except weird indel counter
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
#final table output
with open("loci.tab", "w") as outfile:
	for lc, lcv in sorted(loci.items(), key=operator.itemgetter(1)):
		print >> outfile, lc, "\t", lcv
print "Done"