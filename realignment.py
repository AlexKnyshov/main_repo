#use mafft over a loop for all files in "reduced" directory
import os
import Bio
import sys
from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO
from Bio import AlignIO

filename=sys.argv[1]
files = os.listdir(filename)
if not os.path.exists ("./realigned"):
	os.makedirs("./realigned")
for file in files:
	mafft_cline = MafftCommandline(input="./reduced/"+file)
	outputfile=open ("./realigned/"+file.split(".")[0]+".realigned.fas", "w")
	stdout, stderr = mafft_cline()
	align = AlignIO.read(StringIO(stdout), "fasta")
	outputfile.write(stdout)
	outputfile.close()
