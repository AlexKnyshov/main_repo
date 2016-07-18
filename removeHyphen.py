from Bio import SeqIO
import glob
import sys
import os
if len(sys.argv) == 3:
	if sys.argv[2] == "-1":
		infile = sys.argv[1]
		files = [infile]
	elif sys.argv[2] == "-m":
		files = glob.glob(sys.argv[1]+"/*")
else:
	print "removeHyphen.py script for dealigning sequence files in FASTA format"
	print "existing files will be replaced with new content!"
	print "-----------folder input------------"
	print "FORMAT: python removeHyphen.py [inputfolder] -m"
	print "EXAMPLE: python removeHyphen.py ./fasta/ -m"
	print "------------file input-------------"
	print "FORMAT: python removeHyphen.py [inputfile] -1"
	print "EXAMPLE: python removeHyphen.py ./test.fas -1"
	sys.exit()

for f in files:
	inputf = SeqIO.parse(f, "fasta")
	fnew = f.split("/")
	print "processing file", f
	d = {}
	for seq in inputf:
		line = str(seq.seq).replace("-", "").upper().replace("N", "")
		seq.seq = line
		d[seq.id] = seq.seq
	outupf = open(f, "w")
	for key, value in d.items():
		print >> outupf, ">"+key, "\n", value
	outupf.close()
print "done"