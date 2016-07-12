from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
import re
import sys
import glob
import os
def most_common(lst):
    return max(set(lst), key=lst.count)


if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	opt = sys.argv[2]
else:
	print "FORMAT: python translator.py [folder with fasta] [option: -t, -orf]"
	print "EXAMPLE: python translator.py ./fasta -t"
	sys.exit()

badcount = []
count = 0
print "input folder", inputfolder
files = glob.glob(inputfolder+"/*.fas")
if not os.path.exists ("./translated") and opt == "-t":
	os.makedirs("./translated")
elif opt == "-orf":
	outfile = open("frames.tab", "w")
for f in files:
	newframe = True
	framelist = []
	badframe = []
	transseq = ""
	nuclseq = ""
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
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
				if transseq[:-1].count("*") > 1: #or transseq[:-1].count("X") > 1: #check frame 1
					#print "Nooo"
					nuclseq = nuclseq[1:]+Seq("N", Gapped(IUPAC.ambiguous_dna))
					transseq = nuclseq.translate() #frame2
					#if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 2
					if transseq[:-1].count("*") > 1: #or transseq[:-1].count("X") > 1:
						nuclseq = nuclseq[1:]+Seq("N", Gapped(IUPAC.ambiguous_dna))
						transseq = nuclseq.translate() #frame3
						#if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 3
						if transseq[:-1].count("*") > 1: #or transseq[:-1].count("X") > 1:
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
				if transseq[:-1].count("*") > 1: #or transseq[:-1].count("X") > 1:
					#print "bad seq, reset"
					newframe = True
					goodcount = 0
				else:
					newframe = False
			framelist.append(goodcount)
		# 	print >> debug, seq.id, goodcount
		# else:
		# 	print >> debug, seq.id, "NNN seq"
	goodcount = most_common(framelist)
	if float(len(badframe))/len(alignment) > 0.5:
		#print "bad locus", len(alignment), len(badframe)
		# print >> debug, "BAD LOCUS"
		badcount.append(f)
		goodlocus = False
	else:
		goodlocus = True
		count +=1
	prog = "working on file "+str(f)+": status "+str(goodlocus)+", frame "+str(goodcount)
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
	if goodlocus == True:
		goodlocus = False
		if opt == "-t":
			fnew = f.split("/")
			fn = fnew[len(fnew)-1]
			fn2 = "./translated/"+fn.split(".")[0]+".fas"
			outfile = open(fn2, "w")
			for seq in alignment:
				t1 = str(seq.seq).replace("-", "N")
				t2 = t1.replace("?", "N")
				seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
				if goodcount == 0:
					transseq = seq.seq.translate()
				if goodcount == 1:
					transseq = seq.seq[1:].translate()
				if goodcount == 2:
					transseq = seq.seq[2:].translate()
				unambig_translate = str(transseq).replace("B", "X")
				unambig_translate = unambig_translate.replace("Z", "X")
				unambig_translate = unambig_translate.replace("J", "X")
				unambig_translate = unambig_translate.replace("U", "X")
				unambig_translate = unambig_translate.replace("O", "X")
				print >> outfile, ">"+seq.id, "\n", unambig_translate
			outfile.close()
		elif opt == "-orf":
			print >> outfile, f, "\t", goodcount
if opt == "-orf":
	outfile.close()
print ""
print "good loci:", count
print "bad loci:", badcount
for x in badcount:
	print x
print "done"