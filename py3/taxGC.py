import sys
from Bio import AlignIO
from Bio.SeqUtils import GC

if len(sys.argv) == 2:
	infilename = sys.argv[1]
else:
	print ("FORMAT: python taxGC.py file")
	print ("EXAMPLE: taxGC.py fasta.fas")
	sys.exit()

maindict = {}
fhandle = open(infilename, "r")
for seq in AlignIO.read(fhandle, "fasta"):
	datalist = []
	datalist.append(str(round(GC(str(seq.seq)),2)))
	for pos in range(3):
		datalist.append(str(round(GC(str(seq.seq)[pos::3]),2)))
	maindict[seq.id] = datalist

outhandle = open("GC_per_taxon.csv", "w")
print ("taxon,total,pos1,pos2,pos3", file=outhandle)
for taxon in maindict:
	print (taxon+","+",".join(maindict[taxon]), file=outhandle)
outhandle.close()