from Bio import SeqIO
import sys

infile = sys.argv[1]
inrange = int(sys.argv[2])

with open(infile) as inhandle:
	with open(infile+".edited", "w") as outhandle:
		seqs = SeqIO.parse(inhandle, "fasta")
		counter = 1
		for seq in seqs:
			if counter == inrange:
				print >> outhandle, ">"+seq.id
				print >> outhandle, seq.seq
				break
			counter += 1