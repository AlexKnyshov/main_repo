from Bio import SeqIO
import glob
import os
import shutil
import sys

fold = sys.argv[1]
loci = {}

seq1 = SeqIO.parse(fold, "fasta")
for seq in seq1:
	loci[int(seq.id.split(".")[0][1:])] = int(seq.id.split(".")[-1])

out = open("I13431_nMapped.txt", "w")
for locus, val in sorted(loci.items()):
	for v in range(val):
		print >> out, str(locus)+"\t"+str(v+1)+"\t"+"100"
out.close()

