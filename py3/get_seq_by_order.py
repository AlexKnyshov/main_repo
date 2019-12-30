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
			print (">"+seq.id, file=outhandle)
			print (seq.seq, file=outhandle)
			outhandle.close()
			break
		counter += 1