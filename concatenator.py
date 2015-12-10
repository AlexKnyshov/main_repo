from __future__ import print_function
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Nexus import Nexus
import sys
import glob
import os

inputfolder = sys.argv[1] #folder with phylip
outputfolder = sys.argv[2] #folder for nex files
filterfile = sys.argv[3] #filterlist with names of loci
filesconv = glob.glob(inputfolder+"/*.phylip")
#extlen = len(inputext)

print ("parsing the filter")
filterlist = []
filefilt = open(filterfile, "r")
for f in filefilt:
	filterlist.append(f.strip())
filefilt.close()
print ("done")
print ("converting to nexus")
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
for f in filesconv:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	if fn in filterlist:
		print (fn, end='')
		print('\r' * len(fn), end='')
		alignment = AlignIO.read(inputfolder+"/"+fn, "phylip-relaxed", alphabet=Gapped(IUPAC.ambiguous_dna))
		g = open(outputfolder+"/"+fn[:(len(fn)-6)]+"nex", "w")
		g.write (alignment.format("nexus"))
		g.close()
print('\r' * len(fn), end='')
print("done              ")
print ("concatenation")
filesconcat = glob.glob(outputfolder+"/*.nex")
nexi = []
for f in filesconcat:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	nexi.append((fn, Nexus.Nexus(f)))

 
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('COMBINED.nex', 'w'))
