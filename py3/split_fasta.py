from Bio import SeqIO
import glob
import os
import shutil
import sys
if len(sys.argv) == 2:
    fname = sys.argv[1]
else:
    print ("FORMAT: python split_fasta.py [file name]")
    print ("EXAMPLE: python split_fasta.py ./fasta.fas")
    sys.exit()

fhandle = open(fname, "rU")
seqs = SeqIO.parse(fhandle, "fasta")
for seq in seqs:
    #seq.id = seq.id.split("-")[1]
    print ("processing", seq.id)
    fnew = open(seq.id+".fas", "a")
    SeqIO.write(seq, fnew, "fasta")
    fnew.close()
print ("done")