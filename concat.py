from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
#from Bio.Alphabet import generic_protein
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
import csv
import sys
import glob
import os

def most_common(lst):
    return max(set(lst), key=lst.count)

if len(sys.argv) >= 5:
	inputfolder = sys.argv[1]
	partnum = sys.argv[2]
	phyliptype = sys.argv[3]
	pf2opt = sys.argv[4]
	if partnum == "-orf":
		framefile = sys.argv[5]
else:
	print "FORMAT: python concat.py [folder with fasta] [split to codon positions: -3 (yes), -1 (no), -prot, -orf] [phylip type: -i (interleaved), -s (sequential)] [partition_finder output: -pf2y, -pf2n] ([exclusion list])"
	print "EXAMPLE: python concat.py ./fasta -1 -i -pf2n"
	print "EXAMPLE: python concat.py ./fasta -1 -s -pf2y list.lst"
	print "output is written to COMBINED.phy, partitions are written to partitions.prt"
	sys.exit()

if partnum == "-orf":
	print "reading framefile..."
	loci = {}
	with open(framefile, "rb") as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			#print row[0], row[1]
			loci[row[0].strip()] = int(row[1])

print "creating a list of taxa..."
d = {}
files = glob.glob(inputfolder+"/*.fas")
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	if partnum == "-prot":
		alignment = AlignIO.read(f, "fasta", alphabet = Gapped(IUPAC.protein, '-'))
	else:
		alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	length = alignment.get_alignment_length()
	for seq in alignment:
		if seq.id in d:
			d[seq.id].append(fn)
		else:
			d[seq.id] = []
			d[seq.id].append(fn)
print len(d), "taxa found in fasta alignments"

print "concatenating..."
oks = 0
skips = 0
if partnum == "-prot":
	align = MultipleSeqAlignment([], alphabet = Gapped(IUPAC.protein, '-')) # new ali
else:
	align = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna)) # new ali
for dx in d.keys():
	if partnum == "-prot":
		align.append(SeqRecord(Seq("", alphabet = Gapped(IUPAC.protein, '-'), id=dx))) #add taxa
	else:
		align.append(SeqRecord(Seq("", Gapped(IUPAC.ambiguous_dna)), id=dx)) #add taxa
#print align
counter = 0
outputfile = open("partitions.prt", "w")
if pf2opt == "-pf2y":
	pf2cfg = open("partition_finder.cfg", "w")
 	print >> pf2cfg, "alignment = COMBINED.phy;"
 	print >> pf2cfg, "branchlengths = linked;"
 	print >> pf2cfg, "models = all;"
 	print >> pf2cfg, "model_selection = BIC;"
 	print >> pf2cfg, "[data_blocks]"
#debug = open("debug.log", "w")
badlocus = []
for f in files:
	if partnum == "-prot":
		temp = MultipleSeqAlignment([], alphabet = Gapped(IUPAC.protein, '-')) #temp ali
	else:
		temp = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna)) #temp ali
 	fnew = f.split("/")
 	fn = fnew[len(fnew)-1]
 	if partnum == "-prot":
 		alignment = AlignIO.read(f, "fasta", alphabet = Gapped(IUPAC.protein, '-')) #read original
 	else:
 		alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna)) #read original
 	length = alignment.get_alignment_length()
 	if partnum == "-orf":
	 	remainder = length % 3
		if remainder > 0:
			length = length + (3-remainder)
	missed = list(d)
	goodcount = 0
	
	#print >> debug, "---------------------------------------------------------"
	#print >> debug, f
	if partnum == "-orf":
		if f in loci.keys():
			frame = loci[f]
		else:
			frame = 0
 	for seq in alignment:
 		if partnum == "-orf":
	 		seq.seq = seq.seq[frame:]+Seq("N"*((length - len(seq.seq))+frame), Gapped(IUPAC.ambiguous_dna))
	 	else:
	 		seq.seq = seq.seq
 		#print len(seq.seq)
 		#print seq.id
 		temp.append(seq) #add original to temp
 		if seq.id in missed:
 			missed.remove(seq.id)
 	#print length
 	for m in missed:
 		if partnum == "-prot":
 			temp.append(SeqRecord(Seq("X"*length, alphabet = Gapped(IUPAC.protein, '-'), id=m))) #add dummies
 		else:
 			temp.append(SeqRecord(Seq("?"*(length), Gapped(IUPAC.ambiguous_dna)), id=m)) #add dummies
 	counter = 0
 	if partnum == "-prot":
 		temp2 = MultipleSeqAlignment([], alphabet = Gapped(IUPAC.protein, '-'))
 	else:
 		temp2 = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna))
 	for aliseq in align:
 		for tempseq in temp:
 			if aliseq.id == tempseq.id:
 				temp2.append(aliseq + tempseq)
 	start = align.get_alignment_length()+1
 	end = align.get_alignment_length()+length
 	prog = "working on partition "+str(fn)+": starts "+str(start)+", ends "+str(end)
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
 	align = temp2
 	counter += align.get_alignment_length()
 	if pf2opt == "-pf2y":
 		if partnum == "-3":
			print >> pf2cfg, fn[:-4]+"_1 = "+str(start)+" - "+str(end)+"\\3;"
			print >> pf2cfg, fn[:-4]+"_2 = "+str(start+1)+" - "+str(end)+"\\3;"
			print >> pf2cfg, fn[:-4]+"_3 = "+str(start+2)+" - "+str(end)+"\\3;"
		elif partnum == "-1":
			print >> pf2cfg, fn[:-4]+" = "+str(start)+" - "+str(end)+";"
		elif partnum == "-prot":
			print >> pf2cfg, fn[:-4]+" = "+str(start)+" - "+str(end)+";"
		elif partnum == "-orf":
			print >> pf2cfg, fn[:-4]+" = "+str(start)+" - "+str(end)+";"
 	if partnum == "-3":
		print >> outputfile, "DNA, "+fn+"-1 = "+str(start)+" - "+str(end)+"\\3"
		print >> outputfile, "DNA, "+fn+"-2 = "+str(start+1)+" - "+str(end)+"\\3"
		print >> outputfile, "DNA, "+fn+"-3 = "+str(start+2)+" - "+str(end)+"\\3"
	elif partnum == "-1":
	 	print >> outputfile, "DNA, "+fn+" = "+str(start)+" - "+str(end)
	elif partnum == "-prot":
		print >> outputfile, "WAG, "+fn+" = "+str(start)+" - "+str(end)
	if partnum == "-orf":
		if f not in loci.keys():
			print >> outputfile, "DNA, "+fn+" = "+str(start)+" - "+str(end)
		else:
			print >> outputfile, "DNA, "+fn+"-1-2 = "+str(start)+" - "+str(end)+"\\3, "+str(start+1)+" - "+str(end)+"\\3"
			print >> outputfile, "DNA, "+fn+"-3 = "+str(start+2)+" - "+str(end)+"\\3"
# print >> debug, "---------------------------------------------------------"
# print >> debug, "total bad loci:", len(badlocus)
# for x in badlocus:
# 	print >> debug, x
if pf2opt == "-pf2y":
	print >> pf2cfg, "[schemes]"
	print >> pf2cfg, "search = greedy;"
	pf2cfg.close()
outputfile.close()
print "\ndone"

print "writing..."
if phyliptype == "-i":
	AlignIO.write(align, "COMBINED.phy", "phylip-relaxed")
elif phyliptype == "-s":
	outf = open("COMBINED.phy", "w")
	print >> outf, str(len(align))+" "+str(align.get_alignment_length())
	for seq in align:
		print >> outf, str(seq.id)+" "+str(seq.seq)
	outf.close()
print "done"
