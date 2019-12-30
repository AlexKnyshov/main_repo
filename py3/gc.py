import os
import sys
import glob
import numpy
import csv
from Bio import AlignIO
from Bio.SeqUtils import GC
#filepath input
if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*")
	framefile = sys.argv[2]
else:
	print ("FORMAT: python gc.py [folder with files] [tab file with frames]")
	print ("EXAMPLE: python gc.py ./folder frames.tab")
	sys.exit()
if len(files) == 0:
	print ("no files in the directory")

#starting to process files
loci = {}
progbarc = 0

with open(framefile, "rb") as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		loci[row[0].strip()] = int(row[1])


for f in files:
	infile = open(f, "r")
	alignment = AlignIO.read(infile, "fasta")
	seqlist = []
	if f in loci:
		print ("file", f, "frame", loci[f])
		for seq in alignment:
			stringseq = str(seq.seq)[2+loci[f]::3]
			seq.seq = stringseq
			seqlist.append(GC(seq.seq))
		loci[f] = numpy.mean(seqlist)
		print ("gc content:", loci[f])
	else:
		print ("file", f, "not in the tab file, skipping...")

outf = "locigc3.tab"
with open(outf, "w") as outfile:
	for lc, lcv in loci.items():
		print (lc, "\t", lcv, file=outfile)
print ("output is written to", outf)
print ("done")