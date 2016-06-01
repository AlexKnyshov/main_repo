from Bio import SeqIO

import sys
import glob
import os
import shutil

if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*.fas")
	exclusion_file = sys.argv[2]
else:
	print "FORMAT: python removeTaxa.py [folder with fasta] ([exclusion list])"
	print "EXAMPLE: python removeTaxa.py ./fasta"
	print "EXAMPLE: python removeTaxa.py ./fasta list.lst"
	sys.exit()

print "reading exclusion list..."
exclusion_list = []
exfile = open(exclusion_file, "r")
for line in exfile:
	l = line.strip()
	exclusion_list.append(l)
exfile.close()
print "read", len(exclusion_list), "records"

print "creating an output folder..."
if not os.path.exists ("./reduced/"):
    os.makedirs("./reduced") #creating folder if necessary
else:
    shutil.rmtree("./reduced/") #removing old files
    os.makedirs("./reduced")

print "parsing the files..."
for f in files:
	fn = f.split("/")[-1]
	prog = "working on file "+fn
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
	outputfile=open("./reduced/"+fn, "w")
	for seq in SeqIO.parse(f, "fasta"):
		if seq.id not in exclusion_list:
			print >> outputfile, ">"+seq.id, "\n", seq.seq
	outputfile.close()
print "\ndone"
