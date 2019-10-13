from Bio import SeqIO
import sys


if len(sys.argv) >= 2:
	infname = sys.argv[1]

else:
	print "FORMAT: argument1 = probe file in fasta format"
	print "EXAMPLE: ./fasta.fas"
	print "output is written to [fasta.fas]_mod.fas"
	sys.exit()

with open(infname) as infhandle:
	with open(infname+"_mod.fas", "w") as outh:
		seqs = SeqIO.parse(infhandle, "fasta")
		for seq in seqs:
			print >> outh, ">"+seq.id+"_p1"
			print >> outh, seq.seq