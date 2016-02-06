from Bio import SeqIO
import sys
import glob
import os
if len(sys.argv) == 3:
	if sys.argv[1] == "-1":
		infile = sys.argv[2]
		files = [infile]
	elif sys.argv[1] == "-m":
		files = glob.glob(sys.argv[2]+"/*")
else:
	print "removeR.py script for renaming reversed sequences in FASTA files after mafft alignment"
	print "existing files will be replaced with new content!"
	print "-----------folder input------------"
	print "FORMAT: python removeR.py -m [inputfolder]"
	print "EXAMPLE: python removeR.py -m ./fasta/"
	print "------------file input-------------"
	print "FORMAT: python removeR.py -1 [inputfile]"
	print "EXAMPLE: python removeR.py -1 ./test.fas"
	sys.exit()
for f in files:
	print "processing file", f
	inhandle = open(f, "rU")
	align = SeqIO.parse(inhandle, "fasta")
	d = {}
	for seq in align:
		if "_R_" in seq.id:
			d[seq.id[3:]] = seq.seq
		else:
			d[seq.id] = seq.seq
	inhandle.close()
	outhandle = open(f, "w")
	for key, value in d.items():
		print >> outhandle, ">"+key, "\n", value
	outhandle.close()
