from Bio import SeqIO
import sys
import glob
import os
if len(sys.argv) == 4:
	f = sys.argv[1]
	ext = sys.argv[2]
	files = glob.glob(f+"/*"+ext)
	fmt = sys.argv[3]
elif len(sys.argv) == 3:
	f = sys.argv[1]
	files = [f]
	fmt = sys.argv[2]
else:
	print "FORMAT (single file mode): python taxalist.py [fasta file] [format]"
	print "EXAMPLE: python taxalist.py fasta.fas fasta"
	print "FORMAT (folder mode): python taxalist.py [folder] [extension] [format]"
	print "EXAMPLE: python taxalist.py fastafolder/ .fas fasta"
	sys.exit()

d = {}
for infile in files:
	input_handle = open(infile, "rU")
	alignments = SeqIO.parse(input_handle, fmt)
	for seq in alignments:
		if seq.id in d:
			d[seq.id].append(infile.split("/")[-1])
		else:
			d[seq.id] = []
			d[seq.id].append(infile.split("/")[-1])
	input_handle.close()

for key, value in sorted(d.items()):
	print key+"\t"+str(len(value))+"\t"+str(value)