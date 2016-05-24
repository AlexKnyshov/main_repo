from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment

import sys
import glob
import os

if len(sys.argv) > 1:
	inputfolder = sys.argv[1]
else:
	print "FORMAT: python concat.py [folder with fasta] [split to codon positions: -3 (yes), -1 (no)] [partition_finder output: -pf2y, -pf2n] [phylip type: -i (interleaved), -s (sequential)]"
	print "EXAMPLE: python concat.py ./fasta -1 -i -pf2n"
	print "EXAMPLE: python concat.py ./fasta -1 -s -pf2y"
	print "output is written to COMBINED.phy, partitions are written to partitions.prt"
	sys.exit()

print "input folder", inputfolder
files = glob.glob(inputfolder+"/*.fas")
for f in files:
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	for seq in alignment:
		if seq.seq[0] != "-" and seq.seq[0] != "?" and seq.seq[0] != "N":
			print seq.seq
			break



print "done"
