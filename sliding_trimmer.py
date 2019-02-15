import os
import sys
import glob
from Bio import AlignIO

infile = open(sys.argv[1], "r") #input alignment
window = int(sys.argv[2]) #detection window -- tested 20
step = int(sys.argv[3]) #stepsize -- tested 20
thresh = int(sys.argv[4]) #how many bases in window shoudl be bad - probs need half bad? -- tested 10

inputalignment = AlignIO.read(infile, "fasta")
infile.close()

outname = sys.argv[1].split("/")[-1]

length = inputalignment.get_alignment_length()
outdict = {}
outdict2 = {}
for tx in inputalignment:
	outdict[tx.id] = []
	outdict2[tx.id] = []
	#print outdict[tx.id][0]

#build profile at 40% threshold
for x in range(length):
	col = inputalignment[:,x].lower()
	a = col.count('a')
	t = col.count('t')
	g = col.count('g')
	c = col.count('c')
	tot = a+t+g+c
	if tot > 0:
		ok_bases = set()
		ok_bases.add("-")
		if float(a) / tot > 0.4:
			ok_bases.add("a")
		if float(t) / tot > 0.4:
			ok_bases.add("t")
		if float(g) / tot > 0.4:
			ok_bases.add("g")
		if float(c) / tot > 0.4:
			ok_bases.add("c")
		#print ok_bases
		for seq in range(len(col)):
			#print inputalignment[seq].id
			if col[seq] not in ok_bases:
				outdict[inputalignment[seq].id].append("n")
			else:
				outdict[inputalignment[seq].id].append(col[seq])
			outdict2[inputalignment[seq].id].append(col[seq])
	else:
		outdict[inputalignment[seq].id].append(col[seq])
		outdict2[inputalignment[seq].id].append(col[seq])

# outf = open("basemask"+outname, "w")
# for item, val in outdict.items():
# 	print >> outf, ">"+item
# 	print >> outf, "".join(val)
# outf.close()

#per seq masking
for seqid, seqval in outdict.items():
	i = 0
	while i < length:
		if i+window<=length:
			end = i+window
		else:
			end = length
		if "".join(seqval[i:end]).count("n") > thresh:
			for pos in range(i,end):
				outdict2[seqid][pos] = "n"
		i += step


outf2 = open("masked"+outname, "w")
for item1, val1 in outdict2.items():
	print >> outf2, ">"+item1
	print >> outf2, "".join(val1)
outf2.close()
