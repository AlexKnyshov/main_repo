from Bio import SeqIO
import sys
import glob
import os
import operator
if len(sys.argv) >= 3:
	f = sys.argv[1]
	ext = sys.argv[2]
	files = glob.glob(f+"/*"+ext)
	print "multiple file processing"
	if len(sys.argv) == 4:
		wr = True
	else:
		wr = False
elif len(sys.argv) == 2:
	f = sys.argv[1]
	files = [f]
	print "single file processing"
else:
	print "error"
	sys.exit()
#inputformat = sys.argv[2]
#query = sys.argv[2]
#partnum = len(files)
partnum = 3
d = {}
# for x in files:
# 	d["{0}".format(x)]=[]
#print d
for infile in files:
	input_handle = open(infile, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	#print "read", infile
	for seq in alignments:
		#print seq.id
		if seq.id in d:
			d[seq.id].append(infile.split("/")[-1])
		else:
			d[seq.id] = []
			d[seq.id].append(infile.split("/")[-1])
# taxalist = []
result = {}
prog_c = 0
total_c = len(d)
for key, value in sorted(d.items()):
	c = set(value)
	# if len(c) != len(value):
	# 	print "warning", key, len(value), value
	# else:
	# 	print key, len(value), value
	result[key] = len(value)
	prog_c += 1
	if wr:
		print round(float(prog_c) / total_c *100, 2), "%"
		if len(value)<partnum:
		#	print "Warning, taxon", key, "has only", len(value), "partitions: ", value
		#	rm = raw_input("Would you like to remove it?")
		#	if rm == "y" or rm == "Y":
			del d[key]
			print key, "removed"
	else:
		print key, len(value), value
if wr:
	for infile in files:
		input_handle = open(infile, "rU")
		outhandle = open(infile+".fas", "w")
		alignments = SeqIO.parse(input_handle, "fasta")
		for seq in alignments:
			if seq.id in d:
				print >> outhandle, ">"+seq.id
				print >> outhandle, seq.seq
		input_handle.close()
		outhandle.close()
# for key, value in sorted(result.items(), key=operator.itemgetter(0)):
# 	print key, value
# 	for x in value:
# 		taxalist.append()