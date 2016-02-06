from Bio import SeqIO
import sys
import glob
import os
import operator
if len(sys.argv) == 3:
	f = sys.argv[1]
	ext = sys.argv[2]
	files = glob.glob(f+"/*"+ext)
	print "multiple file processing"
elif len(sys.argv) == 2:
	f = sys.argv[1]
	files = [f]
	print "single file processing"
else:
	print "error"
	sys.exit()
#inputformat = sys.argv[2]
#query = sys.argv[2]
partnum = len(files)
d = {}
# for x in files:
# 	d["{0}".format(x)]=[]
#print d
c=0
for infile in files:
	input_handle = open(infile, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	for seq in alignments:
		#print seq.id
		if seq.id in d:
			d[seq.id].append(infile.split("/")[-1])
		else:
			d[seq.id] = []
			d[seq.id].append(infile.split("/")[-1])
	c+=1
# taxalist = []
result = {}
for key, value in sorted(d.items()):
	print key, value
	result[key] = len(value)
# for key, value in sorted(result.items(), key=operator.itemgetter(0)):
# 	print key, value
# 	for x in value:
# 		taxalist.append()