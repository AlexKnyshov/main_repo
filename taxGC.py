import sys
from Bio import AlignIO
from Bio.SeqUtils import GC

if len(sys.argv) == 3:
	infilename = sys.argv[1]
	informat = sys.argv[2]
else:
	print "FORMAT: python taxGC.py file format"
	print "EXAMPLE: taxGC.py fasta.fas fasta"
	sys.exit()

maindict = {}
fhandle = open(infilename, "r")
for seq in AlignIO.read(fhandle, informat):
	datalist = []
	datalist.append(str(round(GC(str(seq.seq).replace("-", "").replace("?", "").upper().replace("N", "")),2)))
	for pos in range(3):
		datalist.append(str(round(GC(str(seq.seq)[pos::3].replace("-", "").replace("?", "").upper().replace("N", "")),2)))
	maindict[seq.id] = datalist

outhandle = open("GC_per_taxon.csv", "w")
print >> outhandle, "taxon,total,pos1,pos2,pos3"
for taxon in maindict:
	print >> outhandle, taxon+","+",".join(maindict[taxon])
outhandle.close()