from Bio import SeqIO
#from Bio.Alphabet import IUPAC, Gapped
#from Bio.Nexus import Nexus
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
	#print "read OK"
	#count = 0
	d = {}
	for seq in align:
		#print seq.id
		if "_R_" in seq.id:
			#print "found!"
			d[seq.id[3:]] = seq.seq
			#print seq.id, seq.seq
			#seq.name = ""
			#seq.description = ""
		else:
			d[seq.id] = seq.seq
		#count +=1
	inhandle.close()
	#print len(d), count
	outhandle = open(f, "w")
	for key, value in d.items():
		print >> outhandle, ">"+key, "\n", value
	#realign = SeqIO.write(align, outhandle, "fasta")
	outhandle.close()
