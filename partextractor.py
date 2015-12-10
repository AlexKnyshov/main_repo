from Bio import AlignIO
from Bio.Nexus import Nexus
input_handle = open('COMBINED.nex', "rU")
output_handle = open("COMBINED.phy", "w")

print "conversion ....."
alignments = AlignIO.parse(input_handle, "nexus")
AlignIO.write(alignments, output_handle, "phylip-relaxed")
output_handle.close()
input_handle.close()
print "conversion done"

print "generating partitionfinder cfg file"
aln = Nexus.Nexus()
aln.read('COMBINED.nex')
d = {}
for n in aln.charsets:
	start = aln.charsets[n][0]
	end = aln.charsets[n][len(aln.charsets[n])-1]
	d[n] = [start, end]
outputfile = open("partition_finder.cfg", "w")
print >> outputfile, "alignment = COMBINED.phy;"
print >> outputfile, "branchlengths = linked;"
print >> outputfile, "models = all;"
print >> outputfile, "model_selection = BIC;"
print >> outputfile, "[data_blocks]"
for x, y in sorted(d.items()):
	#x(len(fn)-extlen)
	key = str(x)
	print >> outputfile, key[:(len(key)-4)]+"_1 = "+str(y[0])+" - "+str(y[1])+"\\3;"
	print >> outputfile, key[:(len(key)-4)]+"_2 = "+str(y[0]+1)+" - "+str(y[1])+"\\3;"
	print >> outputfile, key[:(len(key)-4)]+"_3 = "+str(y[0]+2)+" - "+str(y[1])+"\\3;"
print >> outputfile, "[schemes]"
print >> outputfile, "search = greedy;"
outputfile.close()
print "done"