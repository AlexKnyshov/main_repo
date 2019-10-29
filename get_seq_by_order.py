from Bio import SeqIO
import sys

infile = sys.argv[1]
inrange = int(sys.argv[2])

with open(infile) as inhandle:
	seqs = SeqIO.parse(inhandle, "fasta")
	counter = 1
	for seq in seqs:
		if counter == inrange:
			outhandle = open(infile+".edited", "w")
			print >> outhandle, ">"+seq.id
			print >> outhandle, seq.seq
			outhandle.close()
			break
		counter += 1