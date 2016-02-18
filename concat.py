from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment

import sys
import glob
import os

if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	partnum = sys.argv[2]
else:
	sys.exit()

files = glob.glob(inputfolder+"/*.fas")

outputfile = open("partitions.prt", "w")

d = {}
print "creating a list of taxa"
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	length = alignment.get_alignment_length()
	for seq in alignment:
		if seq.id in d:
			d[seq.id].append(fn)
		else:
			d[seq.id] = []
			d[seq.id].append(fn)
print "concatenation"
oks = 0
skips = 0
align = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna)) # new ali
for dx in d.keys():
	align.append(SeqRecord(Seq("", Gapped(IUPAC.ambiguous_dna)), id=dx)) #add taxa
#print align 
counter = 0
for f in files:
	temp = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna)) #temp ali
 	fnew = f.split("/")
 	fn = fnew[len(fnew)-1]
 	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna)) #read original
 	length = alignment.get_alignment_length()
	missed = list(d)
 	for seq in alignment:
 		temp.append(seq) #add original to temp
 		missed.remove(seq.id)
 	for m in missed:
 		temp.append(SeqRecord(Seq("?"*length, Gapped(IUPAC.ambiguous_dna)), id=m)) #add dummies
 	counter = 0
 	temp2 = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna))
 	for aliseq in align:
 		for tempseq in temp:
 			if aliseq.id == tempseq.id:
 				temp2.append(aliseq + tempseq)
 	start = align.get_alignment_length()+1
 	end = align.get_alignment_length()+length
 	prog = "partition "+str(fn)+" starts "+str(start)+", ends "+str(end)
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
 	align = temp2
 	counter += align.get_alignment_length()
 	if partnum == "3":
		print >> outputfile, "DNA, "+fn+"-1 = "+str(start)+" - "+str(end)+"\\3"
		print >> outputfile, "DNA, "+fn+"-2 = "+str(start+1)+" - "+str(end)+"\\3"
		print >> outputfile, "DNA, "+fn+"-3 = "+str(start+2)+" - "+str(end)+"\\3"
	elif partnum == "1":
	 	print >> outputfile, "DNA, "+fn+" = "+str(start)+" - "+str(end)
print ""
print "writing"
AlignIO.write(align, "COMBINEX.fas", "phylip-relaxed")
print "done"
