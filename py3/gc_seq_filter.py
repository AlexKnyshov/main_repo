import os
import sys
import glob
import numpy
import shutil
from Bio import AlignIO
from Bio.SeqUtils import GC
#filepath input
if len(sys.argv) == 2:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*")
else:
	print ("FORMAT: python gc_seq_filter.py [folder with files]")
	print ("EXAMPLE: python gc_seq_filter.py ./folder")
	sys.exit()
if len(files) == 0:
	print ("no files in the directory")

print ("creating an output folder...")
if not os.path.exists ("./gc_reduced/"):
    os.makedirs("./gc_reduced") #creating folder if necessary
else:
    shutil.rmtree("./gc_reduced/") #removing old files
    os.makedirs("./gc_reduced")

#starting to process files
outf = open("gc_out.tab", "w")
for f in files:
	fl = f.split("/")[-1]
	infile = open(f, "r")
	alignment = AlignIO.read(infile, "fasta")
	allen = float(alignment.get_alignment_length())
	seqlist = {}
	outliers = []
	for seq in alignment:
		if len(str(seq.seq).replace("-","").replace("N","")) / allen > 0.25:
			seqlist[seq.id] = GC(seq.seq)
		else:
			outliers.append(seq.id)
	avg = numpy.mean(seqlist.values())
	for k, v in seqlist.items():
		if abs(v-avg) > 20:
			print ("locus", fl, "outlier", k, "gc", v, "avg", avg, file=outf)
			outliers.append(k)
	aliout = open("./gc_reduced/"+fl, "w")
	for seq in alignment:
		if seq.id not in outliers:
			print (">"+seq.id, "\n", seq.seq, file=aliout)
	aliout.close()
	infile.close()
outf.close()
print ("done")