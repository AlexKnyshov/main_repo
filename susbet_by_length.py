from Bio import SeqIO
import sys

infilename = sys.argv[1]
passing_length = int(sys.argv[2])
with open(infilename) as inhandle:
	with open(infilename+".subset", "w") as outhandle:
		seqs = SeqIO.parse(inhandle, "fasta")
		for seq in seqs:
			if len(seq.seq) >= passing_length:
				SeqIO.write(seq, outhandle, "fasta")