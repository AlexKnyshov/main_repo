from Bio import SeqIO
import sys
import glob
import os
if len(sys.argv) == 3:
	f = sys.argv[1]
	ext = sys.argv[2]
	files = glob.glob(f+"/*"+ext)
elif len(sys.argv) == 2:
	f = sys.argv[1]
	files = [f]
else:
	print "FORMAT (single file mode): python taxalist.py [fasta file]"
	print "EXAMPLE: python taxalist.py fasta.fas"
	print "FORMAT (single file mode): python taxalist.py [folder] [extension]"
	print "EXAMPLE: python taxalist.py fastafolder/ .fas"
	sys.exit()

d = {}
for infile in files:
	input_handle = open(infile, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	for seq in alignments:
		if seq.id in d:
			d[seq.id].append(infile.split("/")[-1])
		else:
			d[seq.id] = []
			d[seq.id].append(infile.split("/")[-1])
	input_handle.close()

for key, value in sorted(d.items()):
	print key+"\t"+str(len(value))+"\t"+str(value)