from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment

import sys
import os

if len(sys.argv) == 4 or len(sys.argv) == 5:
	infile = sys.argv[1]
	fmt = sys.argv[2]
	tp = sys.argv[3]
	if len(sys.argv) == 5 and sys.argv[4] == "-v":
		rev = True
	else:
		rev = False
else:
	print "FORMAT: python align_editor.py [alignment or single seq file] [format] [type: DNA, Prot] ([-v (retain selected region insted of cutting it out)])"
	print "EXAMPLE: python align_editor.py file.phy phylip-relaxed DNA"
	print "EXAMPLE: python align_editor.py file.phy phylip-relaxed Prot -v"
	sys.exit()

if tp == "DNA":
	alph = Gapped(IUPAC.ambiguous_dna)
elif tp == "Prot":
	alph = Gapped(IUPAC.protein, '-')

alignment = AlignIO.read(infile, fmt, alphabet=alph)

rms = raw_input("Enter start")
rme = raw_input("Enter end")
crop_start = int(rms)
crop_end = int(rme)
if rev:
	alignment = alignment[:, crop_start:crop_end]
else:
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

AlignIO.write(alignment, infile+".edited", fmt)

print "done"