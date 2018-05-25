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
	prop_good_per_frame = []
	stoplist = []
	num_ext_stopslist = []
	for frame in range(0,6):
		counter = 0
		num_stops = 0
		num_ext_stops = 0
		for seq in alignment:
			t1 = str(seq.seq).replace("-", "N")
			t2 = t1.replace("?", "N")
			if frame < 3:
				seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
				#if seq.seq[0:3] != "NNN":
				nuclseq = seq.seq
				remainder = len(nuclseq) % 3
				if remainder > 0:
					nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				nuclseq = nuclseq[frame:]+Seq("N"*frame, Gapped(IUPAC.ambiguous_dna))
				transseq = nuclseq.translate() #frame1
				if transseq[:-1].count("*") > 0: #check number of stops
					num_stops += transseq[:-1].count("*")
					num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
				else:
					counter += 1
					num_stops += transseq[:-1].count("*") # still count #stops to check
					num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			#revcom
			else:
				seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
				nuclseq = seq.seq
				remainder = len(nuclseq) % 3
				if remainder > 0:
					nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				nuclseq = nuclseq[frame:]+Seq("N"*frame, Gapped(IUPAC.ambiguous_dna))
				transseq = nuclseq.translate() #frame4
				if transseq[:-1].count("*") > 0: #check number of stops
					num_stops += transseq[:-1].count("*")
					num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
				else:
					counter += 1
					num_stops += transseq[:-1].count("*") # still count #stops to check
					num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			# #test
			# length = len(nuclseq)
			# start = int(length*0.2)
			# end = length-start
			# nuclseq = nuclseq[start:end]
			# #test
				#print >> debug, seq.id, frame
			#else:
				#print >> debug, seq.id, "NNN seq"
		prop_good_per_frame.append(float(counter)/len(alignment))
		stoplist.append(num_stops)
		num_ext_stopslist.append(num_ext_stops)
	print >> debug, "prop_good_per_frame", prop_good_per_frame, "best proportion", max(prop_good_per_frame), "best frame", prop_good_per_frame.index(max(prop_good_per_frame))
	print >> debug, "stoplist", stoplist, "least stops", min(stoplist), "best frame", stoplist.index(min(stoplist))
	print >> debug, "num_ext_stopslist", num_ext_stopslist, "least stops", min(num_ext_stopslist), "best frame", num_ext_stopslist.index(min(num_ext_stopslist))

	if max(prop_good_per_frame) < cutoff:#modified condition in case several frames are good
		print >> debug, "BAD LOCUS", stoplist
		badloci[f] = stoplist.index(min(stoplist))
		goodlocus = False
		frame = badloci[f]
	else:
		if prop_good_per_frame.index(max(prop_good_per_frame)) != stoplist.index(min(stoplist)):
			print >> debug, "discrepancy btw prop_good_per_frame and stoplist", prop_good_per_frame.index(max(prop_good_per_frame)), stoplist.index(min(stoplist))
			if prop_good_per_frame [stoplist.index(min(stoplist))] > cutoff:
				goodlocus = True
				frame = stoplist.index(min(stoplist))
			elif 0 in stoplist:
				goodlocus = True
				frame = stoplist.index(min(stoplist))
			else:
				print >> debug, "BAD LOCUS", stoplist
				badloci[f] = stoplist.index(min(stoplist))
				goodlocus = False
				frame = badloci[f]
		elif num_ext_stopslist.index(min(num_ext_stopslist)) != stoplist.index(min(stoplist)):
			print >> debug, "discrepancy btw num_ext_stopslist and stoplist", num_ext_stopslist.index(min(num_ext_stopslist)), stoplist.index(min(stoplist))
			if prop_good_per_frame [num_ext_stopslist.index(min(num_ext_stopslist))] > cutoff:
				goodlocus = True
				frame = num_ext_stopslist.index(min(num_ext_stopslist))
			elif 0 in num_ext_stopslist:
				goodlocus = True
				frame = num_ext_stopslist.index(min(num_ext_stopslist))
			else:
				print >> debug, "BAD LOCUS", num_ext_stopslist
				badloci[f] = num_ext_stopslist.index(min(num_ext_stopslist))
				goodlocus = False
				frame = badloci[f]
		else:
			goodlocus = True
			count +=1
			frame = prop_good_per_frame.index(max(prop_good_per_frame))
	prog = "working on file "+str(f)+": status "+str(goodlocus)+", frame "+str(prop_good_per_frame.index(max(prop_good_per_frame)))
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
 	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
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
				if frame < 3:
					seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
					if frame == 1:
						seq.seq = seq.seq[1:]
					if frame == 2:
						seq.seq = seq.seq[2:]
					remainder = len(seq.seq) % 3
					if remainder > 0:
						seq.seq = seq.seq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				else:
					seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
					if frame == 4:
						seq.seq = seq.seq[1:]
					if frame == 5:
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