from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
import csv
import sys
import glob
import os
import shutil

if len(sys.argv) == 4 or len(sys.argv) == 3:
	concat_name = sys.argv[1]
	partname = sys.argv[2]
	if len(sys.argv) == 4:
		par = sys.argv[3]
	else:
		par = "DNA"
else:
	print "FORMAT: python split_phy.py [concatenated file] [partition file] [option: -uniq (append locus name to taxa names), -n (normal mode, keep taxa names as is) ([type: DNA (default), Prot])"
	print "EXAMPLE: python split_phy.py file.phy file.prt"
	print "EXAMPLE: python split_phy.py file.phy file.prt Prot"
	sys.exit()

######
print "creating an output folder..."
if not os.path.exists ("./splitout/"):
    os.makedirs("./splitout") #creating folder if necessary
else:
    shutil.rmtree("./splitout/") #removing old files
    os.makedirs("./splitout")

print "splitting..."
if par == "DNA":
	alph = Gapped(IUPAC.ambiguous_dna)
elif par == "Prot":
	alph = Gapped(IUPAC.protein, '-')
concat_handle = open(concat_name, "rU")
partition_handle = open(partname, "r")
alignment = AlignIO.read(concat_handle, "phylip-relaxed", alphabet = Gapped(IUPAC.ambiguous_dna))
for line in partition_handle:
	l = line.strip().split("=")
	fname = l[0].split(",")[-1].strip() #partition name
	start = int(l[1].split("-")[0])-1
	end = int(l[1].split("-")[1])
	outfile = open("./splitout/"+fname+".fas", "w")
	for record in alignment:
		if set(record.seq[start:end]) != set(['?']):
			print >> outfile, ">"+record.id, "\n",record.seq[start:end]
	outfile.close()
concat_handle.close()
partition_handle.close()