from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment

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
	if len(sys.argv) == 6:
		exclusion_file = sys.argv[5]
else:
	print "FORMAT: python concat.py [folder with fasta] [split to codon positions: -3 (yes), -1 (no), -prot, -orf] [phylip type: -i (interleaved), -s (sequential)] [partition_finder output: -pf2y, -pf2n] ([exclusion list])"
	print "EXAMPLE: python concat.py ./fasta -1 -i -pf2n"
	print "EXAMPLE: python concat.py ./fasta -1 -s -pf2y list.lst"
	print "output is written to COMBINED.phy, partitions are written to partitions.prt"
	sys.exit()

exclusion_list = []
if len(sys.argv) == 6:
	print "reading exclusion list..."
	exfile = open(exclusion_file, "r")
	for line in exfile:
		l = line.strip()
		exclusion_list.append(l)
	exfile.close()
	print "read", len(exclusion_list), "records"

print "creating a list of taxa..."
d = {}
files = glob.glob(inputfolder+"/*.fas")
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	if partnum == "-prot":
		alignment = AlignIO.read(f, "fasta", alphabet = generic_protein)
	else:
		alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	length = alignment.get_alignment_length()
	for seq in alignment:
		if seq.id in d and seq.id not in exclusion_list:
			d[seq.id].append(fn)
		elif seq.id not in exclusion_list:
			d[seq.id] = []
			d[seq.id].append(fn)
print len(d), "taxa found in fasta alignments"

print "concatenating..."
oks = 0
skips = 0
if partnum == "-prot":
	align = MultipleSeqAlignment([], alphabet = generic_protein) # new ali
else:
	align = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna)) # new ali
for dx in d.keys():
	if partnum == "-prot":
		align.append(SeqRecord(Seq("", alphabet = generic_protein), id=dx)) #add taxa
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
for f in files:
	if partnum == "-prot":
		temp = MultipleSeqAlignment([], alphabet = generic_protein) #temp ali
	else:
		temp = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna)) #temp ali
 	fnew = f.split("/")
 	fn = fnew[len(fnew)-1]
 	if partnum == "-prot":
 		alignment = AlignIO.read(f, "fasta", alphabet = generic_protein) #read original
 	else:
 		alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna)) #read original
 	length = alignment.get_alignment_length()
	missed = list(d)
	goodcount = 0
	badlocus = []
	if partnum == "-orf":
		count = 0
		newframe = True
		framelist = []
		badframe = []
		for seq in alignment:
			t1 = str(seq.seq).replace("-", "N")
			t2 = t1.replace("?", "N")
			seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
			if seq.seq[0:3] != "NNN":
				nuclseq = seq.seq
				#break
				#test all seq trans
				#check the reading frame
				remainder = len(nuclseq) % 3
				#print remainder
				if remainder > 0:
					nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				length = len(nuclseq)
				#print length
				if newframe == True: #new search
					#print "new search"
					transseq = nuclseq.translate() #frame1
					#print "locus", f
					#if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 1
					if transseq[:-1].count("*") > 1 or transseq[:-1].count("X") > 1: #check frame 1
						#print "Nooo"
						nuclseq = nuclseq[1:]+Seq("N", Gapped(IUPAC.ambiguous_dna))
						transseq = nuclseq.translate() #frame2
						#if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 2
						if transseq[:-1].count("*") > 1 or transseq[:-1].count("X") > 1:
							nuclseq = nuclseq[1:]+Seq("N", Gapped(IUPAC.ambiguous_dna))
							transseq = nuclseq.translate() #frame3
							#if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 3
							if transseq[:-1].count("*") > 1 or transseq[:-1].count("X") > 1:
								#print "bad seq"
								#bad frame for this seq
								newframe = True
								badframe.append(1)
								goodcount = 0
							else:
								goodcount = 2
								newframe = False
						else:
							goodcount = 1
							newframe = False
					else:
						goodcount = 0
						newframe = False
				elif newframe == False: #old search
					#print "old search"
					nuclseq = nuclseq[goodcount:]+Seq("N"*goodcount, Gapped(IUPAC.ambiguous_dna))
					transseq = nuclseq.translate()
					#if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame
					if transseq[:-1].count("*") > 1 or transseq[:-1].count("X") > 1:
						#print "bad seq, reset"
						newframe = True
						goodcount = 0
					else:
						newframe = False
				framelist.append(goodcount)
		goodcount = most_common(framelist)
		#print "final frame", goodcount
		if len(alignment) - len(badframe) < 1:
			#print "bad locus", len(alignment), len(badframe)
			badlocus.append(fn)
	#length = alignment.get_alignment_length()
	#remainder = length % 3
	#print ""
	#print "goodcount", goodcount, "length", length, "locus", fn, "remainder", remainder
 	for seq in alignment:
 		if length > len(seq.seq):
 			seq.seq = seq.seq+Seq("N"*(length - len(seq.seq)), Gapped(IUPAC.ambiguous_dna))
 		seq.seq = seq.seq[goodcount:]+Seq("N"*goodcount, Gapped(IUPAC.ambiguous_dna))
 		#print len(seq.seq)
 		#print seq.id
 		temp.append(seq) #add original to temp
 		if seq.id in missed:
 			missed.remove(seq.id)
 	#print length
 	for m in missed:
 		if partnum == "-prot":
 			temp.append(SeqRecord(Seq("X"*length, alphabet = generic_protein), id=m)) #add dummies
 		else:
 			temp.append(SeqRecord(Seq("?"*(length), Gapped(IUPAC.ambiguous_dna)), id=m)) #add dummies
 	counter = 0
 	if partnum == "-prot":
 		temp2 = MultipleSeqAlignment([], alphabet = generic_protein)
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
		if fn in badlocus:
			print >> outputfile, "DNA, "+fn+" = "+str(start)+" - "+str(end)
		else:
			print >> outputfile, "DNA, "+fn+"-1 = "+str(start)+" - "+str(end)+"\\3, "+str(start+1)+" - "+str(end)+"\\3"
			print >> outputfile, "DNA, "+fn+"-3 = "+str(start+2)+" - "+str(end)+"\\3"
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
