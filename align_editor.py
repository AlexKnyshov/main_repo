from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment

import sys
import os

infile = sys.argv[1]
#crop_start = int(sys.argv[2])
#crop_end = int(sys.argv[3])

alignment = AlignIO.read(infile, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))

rms = raw_input("Enter start")
rme = raw_input("Enter end")
crop_start = int(rms)
crop_end = int(rme)
alignment = alignment[:, :crop_start]+alignment[:, crop_end:]
q = raw_input("More?")
if q == "y" or q == "Y":
	while q == "y" or q == "Y":
		rms = raw_input("Enter start")
		rme = raw_input("Enter end")
		crop_start = int(rms)
		crop_end = int(rme)
		alignment = alignment[:, :crop_start]+alignment[:, crop_end:]
		q = raw_input("More?")

AlignIO.write(alignment, infile+".edited", "fasta")

print "done"