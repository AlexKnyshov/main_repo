from Bio import SeqIO
import re
import glob
import sys
import os
if len(sys.argv) == 3:
	if sys.argv[1] == "-1":
		infile = sys.argv[2]
		files = [infile]
	elif sys.argv[1] == "-m":
		files = glob.glob(sys.argv[2]+"/*")
else:
	print "removeHyphen.py script for dealigning sequence files in FASTA format"
	print "existing files will be replaced with new content!"
	print "-----------folder input------------"
	print "FORMAT: python removeHyphen.py -m [inputfolder]"
	print "EXAMPLE: python removeHyphen.py -m ./fasta/"
	print "------------file input-------------"
	print "FORMAT: python removeHyphen.py -1 [inputfile]"
	print "EXAMPLE: python removeHyphen.py -1 ./test.fas"
	sys.exit()

for f in files:
	inputf = SeqIO.parse(f, "fasta")
	fnew = f.split("/")
	print "processing file", f
	d = {}
	for seq in inputf:
		line = str(seq.seq)
		line = re.sub('[-]', '', line)
		seq.seq = line
		d[seq.id] = seq.seq
	outupf = open(f, "w")
	for key, value in d.items():
		print >> outupf, ">"+key, "\n", value
	outupf.close()
print "done"