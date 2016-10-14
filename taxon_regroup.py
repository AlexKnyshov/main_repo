import sys
import glob
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo

if len(sys.argv) >= 3:
	if sys.argv[1] == "-tree":
		print "tree mode"
		treefile = sys.argv[3]
		tree = Phylo.read(treefile, "newick")
		tree.ladderize()
		order = tree.get_terminals()
		count = 0
		d = {}
		for x in order:
			print x.name
			d[x.name] = count
			count+=1
	elif sys.argv[1] == "-seqlen":
		print "seqlen mode"
else:
	print "FORMAT: python taxon_regroup.py [option: -tree (regroup based on tree topology), -seqlen (regroup based on seqlen)] [folder] ([tree file])"
	print "EXAMPLE: python taxon_regroup.py -tree ./fasta tree.tre"
	print "EXAMPLE: python taxon_regroup.py -seqlen ./fasta"
	sys.exit()

inputfolder = sys.argv[2]

files = glob.glob(inputfolder+"/*.fas")
if not os.path.exists ("./regrouped"):
	os.makedirs("./regrouped")
for f in files:
	fhandle = open(f, "r")
	sortdict = {}
	for seq in SeqIO.parse(fhandle, "fasta"):
		sortdict[seq.id] = seq.seq
	fhandle.close()
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	fn2 = "./regrouped/"+fn.split(".")[0]+".fas"
	fhandle2 = open(fn2, "w")
	if sys.argv[1] == "-seqlen":
		for key in sorted(sortdict, key=lambda value: len(sortdict[value])):# lambda r: len(str(value).replace("-", "").upper().replace("N", ""))):
			print >> fhandle2, ">L"+str(key)
			print >> fhandle2, sortdict[key]
	if sys.argv[1] == "-tree":
		for key in sorted(sortdict, key = lambda r: d[r.id]):# lambda r: len(str(value).replace("-", "").upper().replace("N", ""))):
			print >> fhandle2, ">L"+str(key)
			print >> fhandle2, sortdict[key]
	fhandle2.close()

	# alignment = AlignIO.read(f, "fasta")
	# if sys.argv[1] == "-tree":
	# 	alignment.sort(key = lambda r: d[r.id])
	# elif sys.argv[1] == "-seqlen":
	# 	alignment.sort(key = lambda r: len(str(r.seq).replace("-", "").upper().replace("N", "")))
	# fnew = f.split("/")
	# fn = fnew[len(fnew)-1]
	# fn2 = "./regrouped/"+fn.split(".")[0]+".fas"
	# outfile = open(fn2, "w")
	# AlignIO.write(alignment, outfile, "fasta")
	# outfile.close()