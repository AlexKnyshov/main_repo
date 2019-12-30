import sys
import glob
import os
from Bio import SeqIO
from Bio import AlignIO

if len(sys.argv) == 2:
	inputfolder = sys.argv[1]
else:
	print ("FORMAT: python missing_data.py [folder]")
	print ("EXAMPLE: missing_data.py ./fasta")
	sys.exit()

files = glob.glob(inputfolder+"/*.fas")
taxaset = set()
locilist = []
maindict = {}
for f in files:
	fhandle = open(f, "r")
	sortdict = {}
	form = "fasta"
	for seq in SeqIO.parse(fhandle, form):
		taxaset.add(seq.id)
		#sortdict[seq.id] = len(str(seq.seq).replace("-", "").upper().replace("N", ""))
		sortdict[seq.id] = len(str(seq.seq).replace("-", "").replace("?", "").upper().replace("N", "")) / float(len(str(seq.seq)))
	fhandle.close()
	maindict[f.split("/")[-1]] = sortdict
	locilist.append(f.split("/")[-1])
l = len(maindict)
outhandle = open("missing_data.csv", "w")
print ("taxon,"+",".join(locilist), file=outhandle)
for taxon in taxaset:
	tempstring = taxon
	for locus in locilist:
		if taxon in maindict[locus].keys():
			tempstring += ","+str(maindict[locus][taxon])
		else:
			tempstring += ",0"
	print (tempstring, file=outhandle)
outhandle.close()