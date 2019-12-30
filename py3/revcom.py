from Bio import SeqIO
import os
import sys
if len(sys.argv) == 2:
    fname = sys.argv[1]
else:
    print ("FORMAT: python revcom.py [file name]")
    print ("EXAMPLE: python revcom.py ./fasta.fas")
    sys.exit()

fhandle = open(fname, "rU")
fnew = open(fname+".revcom.fas", "w")
seqs = SeqIO.parse(fhandle, "fasta")
for seq in seqs:
	seq.seq = seq.seq.reverse_complement()
	SeqIO.write(seq, fnew, "fasta")
fnew.close()
fhandle.close()
print ("done")