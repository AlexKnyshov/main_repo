import os
import sys
import glob
import numpy
#import operator
#import socket
from Bio import AlignIO
from Bio.SeqUtils import GC
#filepath input
if len(sys.argv) == 2:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*")
else:
	print "FORMAT: python gc.py [folder with files]"
	print "EXAMPLE: python gc.py ./folder"
	sys.exit()
if len(files) == 0:
	print "no files in the directory"

#starting to process files
loci = {}
progbarc = 0

for f in files:
	prog = "working on file "+f
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
	infile = open(f, "r")
	alignment = AlignIO.read(infile, "fasta")
	seqlist = []
	for seq in alignment:
		seqlist.append(GC(seq.seq))
	loci[f] = numpy.mean(seqlist)

outf = "locigc.tab"
with open(outf, "w") as outfile:
	for lc, lcv in loci.items():
		print >> outfile, lc, "\t", lcv