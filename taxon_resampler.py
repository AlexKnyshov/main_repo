from Bio import SeqIO
import sys
import glob
import os
import operator

print "init..."
folder = sys.argv[1] #folder with mfs
ext = sys.argv[2] #ext of mfs
files = glob.glob(f+"/*"+ext)
infile = sys.argv[3] #old phylip file

print "trying open phylip..."
ref_align = SeqIO.parse(infile, "phylip")

print "parsing mfs..."
for f in files:
	input_handle = open(f, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")

print "done"