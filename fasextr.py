from Bio import SeqIO
import glob
import os
import shutil
import sys
if len(sys.argv) >= 3:
    fasfilename = sys.argv[1]
    query = sys.argv[2:]
else:
    print "FORMAT: python fasextr.py [fasta file] [query ...]"
    print "EXAMPLE: python fasextr.py fasta.fas Locus1"
    sys.exit()

fasfile = open(fasfilename, "rU")
seqs = SeqIO.parse(fasfile, "fasta")
outfile = open("fasextr.fas", "w")
print "search started..."
for seq in seqs:
	if seq.id in query:
		print >> outfile, ">"+seq.id+"\n"+seq.seq
		query.remove(seq.id)
		print seq.id, "found"
	if len(query) == 0:
		break
print "done"