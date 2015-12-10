from Bio import AlignIO
from Bio.Nexus import Nexus
input_handle = open('COMBINED.nex', "rU")
print "generating partition file"
aln = Nexus.Nexus()
aln.read('COMBINED.nex')
d = {}
for n in aln.charsets:
	start = aln.charsets[n][0]
	end = aln.charsets[n][len(aln.charsets[n])-1]
	d[n] = [start, end]
outputfile = open("partitions.prt", "w")
for x, y in sorted(d.items()):
	#x(len(fn)-extlen)
	key = str(x)
	print >> outputfile, "DNA, "+key[:(len(key)-4)]+"-1 = "+str(y[0]+1)+" - "+str(y[1]+1)+"\\3"
	print >> outputfile, "DNA, "+key[:(len(key)-4)]+"-2 = "+str(y[0]+2)+" - "+str(y[1]+1)+"\\3"
	print >> outputfile, "DNA, "+key[:(len(key)-4)]+"-3 = "+str(y[0]+3)+" - "+str(y[1]+1)+"\\3"
outputfile.close()
input_handle.close()
print "done"