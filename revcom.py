from Bio import SeqIO
import os
import sys
if len(sys.argv) == 2:
    fname = sys.argv[1]
else:
    print "FORMAT: python revcom.py [file name]"
    print "EXAMPLE: python revcom.py ./fasta.fas"
    sys.exit()

fhandle = open(fname, "rU")
seq = SeqIO.read(fhandle, "fasta")
seq.seq = seq.seq.reverse_complement()
fnew = open(fname+".revcom.fas", "w")
SeqIO.write(seq, fnew, "fasta")
fnew.close()
fhandle.close()
print "done"