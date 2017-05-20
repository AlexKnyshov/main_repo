from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
import csv
import sys
import glob
import os

concat_name = sys.argv[1]
partname = sys.argv[2]
opt = sys.argv[3]

######
concat_handle = open(concat_name, "rU")
partition_handle = open(partname, "r")
alignment = AlignIO.read(concat_handle, "phylip-relaxed", alphabet = Gapped(IUPAC.protein, '-'))
for line in partition_handle:
	l = line.strip().split("=")
	fname = l[0].split(",")[-1].strip()#, l[1].split(",")
	#print fname
	#outfile = open(fname+".fas", "w")
	i = 1
	for p in l[1].split(","):
		start = int(p.split("-")[0])-1
		end = int(p.split("-")[1])
		outfile = open(fname+"-"+str(i)+".fas", "w")
		for record in alignment:
			#print start, end
			#print record.id, record.seq[start:end]
			if opt == "-uniq":
				print >> outfile, ">"+record.id+"_"+fname+"-"+str(i), "\n",record.seq[start:end]
			else:
				print >> outfile, ">"+record.id, "\n",record.seq[start:end]
		#outfile
		i += 1


concat_handle.close()
partition_handle.close()