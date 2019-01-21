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
###
def alidist(f, framesfunc):
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	distcalc1 = []
	distcalc2 = []
	#print >> debug, "debug framesfunc", framesfunc
	for seq in alignment:
		frame = framesfunc[0]
		t1 = str(seq.seq).replace("-", "N")
		t2 = t1.replace("?", "N")
		if frame < 3:
			nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
			if frame == 1:
				nuclseq = nuclseq[1:]
			if frame == 2:
				nuclseq = nuclseq[2:]
			remainder = len(nuclseq) % 3
			if remainder > 0:
				nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
		else:
			nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
			if frame == 4:
				nuclseq = nuclseq[1:]
			if frame == 5:
				nuclseq = nuclseq[2:]
			remainder = len(nuclseq) % 3
			if remainder > 0:
				nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
		transseq = nuclseq.translate()
		t = str(nuclseq.translate()).upper().replace("*","X")
		transseq = str(transseq)
		distcalc1.append(transseq)

		######
		frame = framesfunc[1]
		if frame < 3:
			nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
			if frame == 1:
				nuclseq = nuclseq[1:]
			if frame == 2:
				nuclseq = nuclseq[2:]
			remainder = len(nuclseq) % 3
			if remainder > 0:
				nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
		else:
			nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
			if frame == 4:
				nuclseq = nuclseq[1:]
			if frame == 5:
				nuclseq = nuclseq[2:]
			remainder = len(nuclseq) % 3
			if remainder > 0:
				nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
		transseq = nuclseq.translate()
		t = str(nuclseq.translate()).upper().replace("*","X")
		transseq = str(transseq)
		distcalc2.append(transseq)
	posx1 = []
	posy1 = 0
	#print >> debug, "debug framesfunc posy1", len(distcalc1), distcalc1, distcalc1[0], distcalc1[0][0]
	for x1 in range(0, len(distcalc1[0])):#pos
		pos1 = []
		for y1 in range(0, len(distcalc1)):#seqs
			pos1.append(distcalc1[y1][x1])
		posx1.append(len(set(pos1)))
	posy1 = sum(posx1)
	posx2 = []
	posy2 = 0
	for x2 in range(0, len(distcalc2[0])):
		pos2 = []
		for y2 in range(0, len(distcalc2)):
			pos2.append(distcalc2[y2][x2])
		posx2.append(len(set(pos2)))
	posy2 = sum(posx2)
	#print >> debug, "debug framesfunc posy1", distcalc1, distcalc1[0]
	if posy1 < posy2:
		#print >> debug, "debug framesfunc posy1", posy1, posy2, posx1, posx2, pos1, pos2
		return framesfunc[0]
	else:
		#print >> debug, "debug framesfunc posy2", posy1, posy2, posx1, posx2, pos1, pos2
		return framesfunc[1]
#######

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
				nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
				if frame == 1:
					nuclseq = nuclseq[1:]
				if frame == 2:
					nuclseq = nuclseq[2:]
				remainder = len(nuclseq) % 3
				if remainder > 0:
					nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
			else:
				nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
				if frame == 4:
					nuclseq = nuclseq[1:]
				if frame == 5:
					nuclseq = nuclseq[2:]
				remainder = len(nuclseq) % 3
				if remainder > 0:
					nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
			transseq = nuclseq.translate()
			if transseq[:-1].count("*") > 0: #check number of stops
				num_stops += transseq[:-1].count("*")
				num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			else:
				counter += 1
				num_stops += transseq[:-1].count("*") # still count #stops to check
				num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			# if frame < 3:
			# 	nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
			# 	#if seq.seq[0:3] != "NNN":
			# 	remainder = len(nuclseq) % 3
			# 	if remainder > 0:
			# 		nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
			# 	nuclseq = nuclseq[frame:]+Seq("N"*frame, Gapped(IUPAC.ambiguous_dna))
			# 	transseq = nuclseq.translate() #frame1
			# 	if transseq[:-1].count("*") > 0: #check number of stops
			# 		num_stops += transseq[:-1].count("*")
			# 		num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			# 	else:
			# 		counter += 1
			# 		num_stops += transseq[:-1].count("*") # still count #stops to check
			# 		num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			# #revcom
			# else:
			# 	nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
			# 	remainder = len(nuclseq) % 3
			# 	if remainder > 0:
			# 		nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
			# 	nuclseq = nuclseq[(frame-3):]+Seq("N"*(frame-3), Gapped(IUPAC.ambiguous_dna))
			# 	transseq = nuclseq.translate() #frame4
			# 	if transseq[:-1].count("*") > 0: #check number of stops
			# 		num_stops += transseq[:-1].count("*")
			# 		num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
			# 	else:
			# 		counter += 1
			# 		num_stops += transseq[:-1].count("*") # still count #stops to check
			# 		num_ext_stops += transseq[int(0.2*len(transseq)):(len(transseq)-int(0.2*len(transseq)))].count("*")
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

	if max(prop_good_per_frame) < cutoff and min(num_ext_stopslist) > 0:#modified condition in case several frames are good
		print >> debug, "BAD LOCUS", stoplist
		badloci[f] = stoplist.index(min(stoplist))
		goodlocus = False
		frame = badloci[f]
	else:
		bestframeslst = [prop_good_per_frame.index(max(prop_good_per_frame)), stoplist.index(min(stoplist)), num_ext_stopslist.index(min(num_ext_stopslist))]
		getindex_1 = [i2 for i2,x in enumerate(prop_good_per_frame) if x == max(prop_good_per_frame)]
		if len(getindex_1) == 2:
			print >> debug, "two identical 1, frames", getindex_1
		getindex_2 = [i2 for i2,x in enumerate(stoplist) if x == min(stoplist)]
		if len(getindex_2) == 2:
			print >> debug, "two identical 2, frames", getindex_2
		getindex_3 = [i2 for i2,x in enumerate(num_ext_stopslist) if x == min(num_ext_stopslist)]
		if len(getindex_3) == 2:
			print >> debug, "two identical 3, frames", getindex_3
		#print >> debug, "getindex", getindex_1, getindex_2, getindex_3
		#option for internal!! correct later
		if min(num_ext_stopslist) == 0 and min(num_ext_stopslist) < min(stoplist) and len(getindex_3) == 1:
			print >> debug, "conflict btw stoplist and ext_stoplist, the latter has a single 0 stop frame", getindex_3[0]
			frame = getindex_3[0]
			goodlocus = True
			count +=1
		else:
			if getindex_1 == getindex_2 and len(getindex_1) == len(getindex_2) == 2:
				print >> debug, "checking average distance, frames", getindex_1
				frame = alidist(f, getindex_1)
				print >> debug, "frame with lowest distance", frame
				goodlocus = True
			else:
				if len(set(bestframeslst)) < 3:
					frame = max(set(bestframeslst), key=bestframeslst.count)
					goodlocus = True
					count +=1
				else:
					badloci[f] = stoplist.index(min(stoplist))
					print >> debug, "BAD LOCUS", bestframeslst
					goodlocus = False
					frame = badloci[f]
		# if prop_good_per_frame.index(max(prop_good_per_frame)) != stoplist.index(min(stoplist)):
		# 	print >> debug, "discrepancy btw prop_good_per_frame and stoplist", prop_good_per_frame.index(max(prop_good_per_frame)), stoplist.index(min(stoplist))
		# 	if prop_good_per_frame [stoplist.index(min(stoplist))] > cutoff:
		# 		goodlocus = True
		# 		frame = stoplist.index(min(stoplist))
		# 	elif 0 in stoplist:
		# 		goodlocus = True
		# 		frame = stoplist.index(min(stoplist))
		# 	else:
		# 		print >> debug, "BAD LOCUS", stoplist
		# 		badloci[f] = stoplist.index(min(stoplist))
		# 		goodlocus = False
		# 		frame = badloci[f]
		# elif num_ext_stopslist.index(min(num_ext_stopslist)) != stoplist.index(min(stoplist)):
		# 	print >> debug, "discrepancy btw num_ext_stopslist and stoplist", num_ext_stopslist.index(min(num_ext_stopslist)), stoplist.index(min(stoplist))
		# 	if prop_good_per_frame [num_ext_stopslist.index(min(num_ext_stopslist))] > cutoff:
		# 		goodlocus = True
		# 		frame = num_ext_stopslist.index(min(num_ext_stopslist))
		# 	elif 0 in num_ext_stopslist:
		# 		goodlocus = True
		# 		frame = num_ext_stopslist.index(min(num_ext_stopslist))
		# 	else:
		# 		print >> debug, "BAD LOCUS", num_ext_stopslist
		# 		badloci[f] = num_ext_stopslist.index(min(num_ext_stopslist))
		# 		goodlocus = False
		# 		frame = badloci[f]
		# else:
		# 	goodlocus = True
		# 	count +=1
		# 	frame = prop_good_per_frame.index(max(prop_good_per_frame))
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
					nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
					if frame == 1:
						nuclseq = nuclseq[1:]
					if frame == 2:
						nuclseq = nuclseq[2:]
					remainder = len(nuclseq) % 3
					if remainder > 0:
						nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				else:
					nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
					if frame == 4:
						nuclseq = nuclseq[1:]
					if frame == 5:
						nuclseq = nuclseq[2:]
					remainder = len(nuclseq) % 3
					if remainder > 0:
						nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				transseq = nuclseq.translate()
				t = str(nuclseq.translate()).upper().replace("*","X")
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


		for frame_num in range(0,6):
			print >> debug, "PREDICTED alignment for frame", frame_num
			for seq in alignment:
				t1 = str(seq.seq).replace("-", "N")
				t2 = t1.replace("?", "N")
				if frame_num < 3:
					nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
					if frame_num == 1:
						nuclseq = nuclseq[1:]
					if frame_num == 2:
						nuclseq = nuclseq[2:]
					remainder = len(nuclseq) % 3
					if remainder > 0:
						nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				else:
					nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
					if frame_num == 4:
						nuclseq = nuclseq[1:]
					if frame_num == 5:
						nuclseq = nuclseq[2:]
					remainder = len(nuclseq) % 3
					if remainder > 0:
						nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
				transseq = nuclseq.translate()
				#t = str(nuclseq.translate()).upper().replace("*","X")
				print >> debug, transseq
				#print >> debug, t
		#alignment.seek(0)
		print >> debug, "---------------------------------------------------------"
if opt == "-orf":
	outfile.close()
print ""
print "good loci:", count
for x, y in badloci.items():
	print x, "predicted frame:", y
print "bad loci:", len(badloci)
print "done"