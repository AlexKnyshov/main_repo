import sys
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Consensus import *
msa = AlignIO.read(sys.argv[1], 'phylip-relaxed')
print "generating bootstrap replicates..."
count = 0
msas = bootstrap(msa, int(sys.argv[2]))
for m in msas:
	count += 1
	print "replicate", count
	print len(m), m.get_alignment_length()
	AlignIO.write(m, sys.argv[1]+".r"+str(count), "phylip-relaxed")
print "done"