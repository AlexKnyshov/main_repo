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


if len(sys.argv) == 4:
	inputfolder = sys.argv[1]
	opt = sys.argv[2]
	cutoff = float(sys.argv[3])
else:
	print "FORMAT: python translator.py [folder with fasta] [option: -t, -orf] [cutoff]"
	print "EXAMPLE: python translator.py ./fasta -t 0.4"
	sys.exit()

badloci = {}
count = 0
print "input folder", inputfolder
files = glob.glob(inputfolder+"/*.fas")
if not os.path.exists ("./translated") and opt == "-t":
	os.makedirs("./translated")
elif opt == "-orf":
	outfile = open("frames.tab", "w")
debug = open("debug.log", "w")
for f in files:
	print >> debug, f
	newframe = True
	transseq = ""
	nuclseq = ""
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	framelist = []
	stoplist = []
	bf2list = []
	for frame in range(0,3):
		counter = 0
		badframe = 0
		bf2 = 0
		for seq in alignment:
			t1 = str(seq.seq).replace("-", "N")
			t2 = t1.replace("?", "N")
			seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
			#if seq.seq[0:3] != "NNN":
			nuclseq = seq.seq
			remainder = len(nuclseq) % 3
			if remainder > 0:
				nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
			nuclseq = nuclseq[frame:]+Seq("N"*frame, Gapped(IUPAC.ambiguous_dna))
			transseq = nuclseq.translate() #frame1
			if transseq[:-1].count("*") > 1: #check number of stops
				badframe += transseq[:-1].count("*")
				bf2 += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			else:
				counter += 1
				badframe += transseq[:-1].count("*") # still count #stops to check
				bf2 += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			# #test
			# length = len(nuclseq)
			# start = int(length*0.2)
			# end = length-start
			# nuclseq = nuclseq[start:end]
			# #test
				#print >> debug, seq.id, frame
			#else:
				#print >> debug, seq.id, "NNN seq"
		framelist.append(float(counter)/len(alignment))
		stoplist.append(badframe)
		bf2list.append(bf2)
	print >> debug, "framelist", framelist, max(framelist), framelist.index(max(framelist))
	print >> debug, "stoplist", stoplist, min(stoplist), stoplist.index(min(stoplist))
	print >> debug, "bf2list", bf2list, min(bf2list), bf2list.index(min(bf2list))

	
	if max(framelist) < cutoff:
		print >> debug, "BAD LOCUS", stoplist
		badloci[f] = stoplist.index(min(stoplist))
		goodlocus = False
		frame = badloci[f]
	else:
		if framelist.index(max(framelist)) != stoplist.index(min(stoplist)):
			print >> debug, "discrepancy btw framelist and stoplist", framelist.index(max(framelist)), stoplist.index(min(stoplist))
			if framelist [stoplist.index(min(stoplist))] > cutoff:
				goodlocus = True
				frame = stoplist.index(min(stoplist))
			else:
				print >> debug, "BAD LOCUS", stoplist
				badloci[f] = stoplist.index(min(stoplist))
				goodlocus = False
				frame = badloci[f]
		elif bf2list.index(min(bf2list)) != stoplist.index(min(stoplist)):
			print >> debug, "discrepancy btw bf2list and stoplist", bf2list.index(min(bf2list)), stoplist.index(min(stoplist))
			if framelist [bf2list.index(min(bf2list))] > cutoff:
				goodlocus = True
				frame = bf2list.index(min(bf2list))
			else:
				print >> debug, "BAD LOCUS", bf2list
				badloci[f] = bf2list.index(min(bf2list))
				goodlocus = False
				frame = badloci[f]
		else:
			goodlocus = True
			count +=1
			frame = framelist.index(max(framelist))
	prog = "working on file "+str(f)+": status "+str(goodlocus)+", frame "+str(framelist.index(max(framelist)))
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
				if frame == 1:
					seq.seq = seq.seq[1:]
				if frame == 2:
					seq.seq = seq.seq[2:]
				remainder = len(seq.seq) % 3
				if remainder > 0:
					seq.seq = seq.seq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				transseq = seq.seq.translate()
				t = str(seq.seq.translate()).upper().replace("*","X")
				print >> debug, transseq
				print >> debug, t
				transseq = str(transseq)
				#if remainder > 0:
				for t in range(len(transseq)-1,-1,-1):
					if transseq[t] != "X":
						transseq = transseq[:t]+"X"+transseq[t+1:]
						break
				#unambig_translate = transseq.replace("B", "X")
				#unambig_translate = unambig_translate.replace("Z", "X")
				#unambig_translate = unambig_translate.replace("J", "X")
				#unambig_translate = unambig_translate.replace("U", "X")
				#unambig_translate = unambig_translate.replace("O", "X")
				print >> outfile, ">"+seq.id, "\n", transseq#unambig_translate
			outfile.close()
		elif opt == "-orf":
			print >> outfile, f, "\t", frame
		print >> debug, "Final frame:", frame
		print >> debug, "---------------------------------------------------------"
	else:
		if opt == "-orf":
			print >> outfile, f, "\t", frame
		print >> debug, "PREDICTED final frame:", frame
		print >> debug, "---------------------------------------------------------"
if opt == "-orf":
	outfile.close()
print ""
print "good loci:", count
for x, y in badloci.items():
	print x, "predicted frame:", y
print "bad loci:", len(badloci)
print "done"