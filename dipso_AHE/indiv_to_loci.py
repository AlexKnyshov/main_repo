from Bio import SeqIO
import glob
import os
import shutil
import sys

fold = sys.argv[1]
files1 = glob.glob(fold+"/*.fasta")
loci = set()

for f in files1:
	seq1 = SeqIO.parse(f, "fasta")
	for seq in seq1:
		loci.add(seq.id)

for locus in loci:
    createfile = open(locus+".fas", "w")
    createfile.close()

for f in files1:
	seq1 = SeqIO.parse(f, "fasta")
	for seq in seq1:
		appendfile = open(seq.id+".fas", "a")
		print >> appendfile, ">"+f.split("/")[-1], "\n", seq.seq
		appendfile.close()
