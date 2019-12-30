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

if len(sys.argv) >= 3:
	inputfolder = sys.argv[1]
	opt = sys.argv[2]
	if opt == "-t" or opt == "-orf":
		cutoff = float(sys.argv[3])
else:
	print ("FORMAT: python translator.py [folder with fasta] [option: -t, -orf, -u] ([cutoff])")
	print ("EXAMPLE: python translator.py ./fasta -t 0.4")
	print ("EXAMPLE: python translator.py ./fasta -u")
	sys.exit()
###
def output_partial(f, frame):
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	distcalc1 = []
	prelist = []
	# for f in alignment:
	# 	prelist.append(SeqRecord(Seq("N"*gap, Gapped(IUPAC.ambiguous_dna, '-')), id=f.id))
	# alignreplacement = MultipleSeqAlignment(prelist)
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
		prelist.append(SeqRecord(transseq, id=seq.id))
		t = str(nuclseq.translate()).upper().replace("*","X")
		transseq = str(transseq)
		distcalc1.append(transseq)
	posx1 = []
	pos_stop = []
	pos_X = []
	#print >> debug, "debug framesfunc posy1", len(distcalc1), distcalc1, distcalc1[0], distcalc1[0][0]
	start = 0
	end = len(distcalc1[0])
	len1 = 0
	region_score = []
	region_coords = []
	for x1 in range(0, len(distcalc1[0])):#pos
		pos1 = []
		for y1 in range(0, len(distcalc1)):#seqs
			pos1.append(distcalc1[y1][x1])
		pos_stop.append(pos1.count("*"))
		if pos_stop[x1] > 0:
			end = x1
			len1 = end - start
			region_score.append(len1)
			region_coords.append((start, end))
			start = x1
			len1 = 1
		else:
			len1 += 1
		pos_X.append(pos1.count("X"))
	end = x1
	len1 = end - start
	region_score.append(len1)
	region_coords.append((start, end))
	# if frame == 1 or frame == 3:
	# 	buffer1 = 1
	# elif frame == 2 or frame == 4:
	# 	buffer1 = 2
	# else:
	# 	buffer1 = 0
	good_coords = region_coords[region_score.index(max(region_score))]
	alignreplacement = MultipleSeqAlignment(prelist)
	print (alignreplacement[:,good_coords[0]+1:good_coords[1]-1], file=debug)
	#print >> debug, (region_coords[region_score.index(max(region_score))][0])*3+buffer1, region_coords[region_score.index(max(region_score))][1]*3
	return alignreplacement[:,good_coords[0]+1:good_coords[1]-1]

def profile_check(f, frame):
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	distcalc1 = []
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
		transseq = str(transseq)
		distcalc1.append(transseq)
	posx1 = []
	pos_stop = []
	pos_X = []
	#print >> debug, "debug framesfunc posy1", len(distcalc1), distcalc1, distcalc1[0], distcalc1[0][0]
	for x1 in range(0, len(distcalc1[0])):#pos
		pos1 = []
		for y1 in range(0, len(distcalc1)):#seqs
			pos1.append(distcalc1[y1][x1])
		posx1.append(len(set(pos1)))
		pos_stop.append(pos1.count("*"))
		pos_X.append(pos1.count("X"))
	while float(sum(pos_stop))/len(alignment) > 0.4:
		print (float(sum(pos_stop))/len(alignment))
		print (sum(pos_stop))
		#try to improve
		if pos_stop[0:int(len(posx1)*0.3)] < pos_stop[(len(posx1)-1):int(len(posx1)*0.7):-1]:
			cumul = pos_stop[0:5]
			for run1 in range(5, len(posx1), 1):
				print (sum(cumul)/5, file=debug)
				if sum(cumul)/5 > 0:
					print (run1, file=debug)
					break
				del cumul[0]
				cumul.append(pos_stop[run1])
		else:
			cumul = pos_stop[len(posx1)-5:len(posx1)]
			for run1 in range(len(posx1)-5, 0, -1):
				#print >> debug, cumul, posx1[run1]
				print (sum(cumul)/5, file=debug)
				if sum(cumul)/5 > 0:
					print (run1, file=debug)
					break
				del cumul[0]
				cumul.append(pos_stop[run1])
		bestframe = frame
		bestframescore = sum(pos_stop)
		for gap in range(1,3):
			prelist = []
			for f in alignment:
				prelist.append(SeqRecord(Seq("N"*gap, Gapped(IUPAC.ambiguous_dna, '-')), id=f.id))
			alignreplacement = MultipleSeqAlignment(prelist)
			alignment = alignment[:,0:run1*3]+alignreplacement+alignment[:,run1*3:]
			for fr in range(0,6):
				distcalc1 = []
				for seq in alignment:
					t1 = str(seq.seq).replace("-", "N")
					t2 = t1.replace("?", "N")
					if fr < 3:
						nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
						if fr == 1:
							nuclseq = nuclseq[1:]
						if fr == 2:
							nuclseq = nuclseq[2:]
						remainder = len(nuclseq) % 3
						if remainder > 0:
							nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
					else:
						nuclseq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna)).reverse_complement()
						if fr == 4:
							nuclseq = nuclseq[1:]
						if fr == 5:
							nuclseq = nuclseq[2:]
						remainder = len(nuclseq) % 3
						if remainder > 0:
							nuclseq = nuclseq+Seq("N"*(3-remainder), Gapped(IUPAC.ambiguous_dna))
					transseq = nuclseq.translate()
					t = str(nuclseq.translate()).upper().replace("*","X")
					transseq = str(transseq)
					distcalc1.append(transseq)
				posx1 = []
				pos_stop = []
				pos_X = []
				#print >> debug, "debug framesfunc posy1", len(distcalc1), distcalc1, distcalc1[0], distcalc1[0][0]
				for x1 in range(0, len(distcalc1[0])):#pos
					pos1 = []
					for y1 in range(0, len(distcalc1)):#seqs
						pos1.append(distcalc1[y1][x1])
					posx1.append(len(set(pos1)))
					pos_stop.append(pos1.count("*"))
					pos_X.append(pos1.count("X"))
				print (run1, file=debug)
				print ("gap", gap, "frame", fr, file=debug)
				print (sum(posx1), posx1, file=debug)
				print (sum(pos_stop), pos_stop, file=debug)
				print (sum(pos_X), pos_X, file=debug)
				print (distcalc1, file=debug)
				if sum(pos_stop) < bestframescore:
					bestframe = fr
					bestframescore = sum(pos_stop)
					bestali = distcalc1
	print ("frame", bestframe, bestframescore, bestali, file=debug)


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
print ("input folder", inputfolder)
files = glob.glob(inputfolder+"/*.fas")
if not os.path.exists ("./translated") and (opt == "-t" or opt == "-u"):
	os.makedirs("./translated")
elif opt == "-orf":
	outfile = open("frames.tab", "w")
debug = open("debug.log", "w")
if opt == "-t" or opt == "-orf":
	for f in files:
		print (f, file=debug)
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
			prop_good_per_frame.append(float(counter)/len(alignment))
			stoplist.append(num_stops)
			num_ext_stopslist.append(num_ext_stops)
		print ("prop_good_per_frame", prop_good_per_frame, "best proportion", max(prop_good_per_frame), "best frame", prop_good_per_frame.index(max(prop_good_per_frame)), file=debug)
		print ("stoplist", stoplist, "least stops", min(stoplist), "best frame", stoplist.index(min(stoplist)), file=debug)
		print ("num_ext_stopslist", num_ext_stopslist, "least stops", min(num_ext_stopslist), "best frame", num_ext_stopslist.index(min(num_ext_stopslist)), file=debug)
		#addopt = False
		if max(prop_good_per_frame) < cutoff and min(num_ext_stopslist) > 0:#modified condition in case several frames are good
			print ("BAD LOCUS", stoplist, file=debug)
			#additional tests:
			print ("TEST", file=debug)
			print ("TEST", file=debug)
			#profile_check(f,stoplist.index(min(stoplist)))
			partial = output_partial(f,stoplist.index(min(stoplist)))
			fnew = f.split("/")
			fn = fnew[len(fnew)-1]
			fn2 = "./translated/"+fn.split(".")[0]+".fas"
			outfile = open(fn2, "w")
			AlignIO.write(partial, outfile, "fasta")
			outfile.close()
			#addopt = True
			badloci[f] = stoplist.index(min(stoplist))
			goodlocus = False
			frame = badloci[f]
			#frame = stoplist.index(min(stoplist))
		else:
			bestframeslst = [prop_good_per_frame.index(max(prop_good_per_frame)), stoplist.index(min(stoplist)), num_ext_stopslist.index(min(num_ext_stopslist))]
			getindex_1 = [i2 for i2,x in enumerate(prop_good_per_frame) if x == max(prop_good_per_frame)]
			if len(getindex_1) == 2:
				print ("two identical 1, frames", getindex_1, file=debug)
			getindex_2 = [i2 for i2,x in enumerate(stoplist) if x == min(stoplist)]
			if len(getindex_2) == 2:
				print ("two identical 2, frames", getindex_2, file=debug)
			getindex_3 = [i2 for i2,x in enumerate(num_ext_stopslist) if x == min(num_ext_stopslist)]
			if len(getindex_3) == 2:
				print ("two identical 3, frames", getindex_3, file=debug)
			#print >> debug, "getindex", getindex_1, getindex_2, getindex_3
			#option for internal!! correct later
			if min(num_ext_stopslist) == 0 and min(num_ext_stopslist) < min(stoplist) and len(getindex_3) == 1:
				print ("conflict btw stoplist and ext_stoplist, the latter has a single 0 stop frame", getindex_3[0], file=debug)
				frame = getindex_3[0]
				goodlocus = True
				count +=1
			else:
				if getindex_1 == getindex_2 and len(getindex_1) == len(getindex_2) == 2:
					print ("checking average distance, frames", getindex_1, file=debug)
					frame = alidist(f, getindex_1)
					print ("frame with lowest distance", frame, file=debug)
					goodlocus = True
				else:
					if len(set(bestframeslst)) < 3:
						frame = max(set(bestframeslst), key=bestframeslst.count)
						goodlocus = True
						count +=1
					else:
						badloci[f] = stoplist.index(min(stoplist))
						print ("BAD LOCUS", bestframeslst, file=debug)
						goodlocus = False
						frame = badloci[f]
		prog = "working on file "+str(f)+": status "+str(goodlocus)+", frame "+str(prop_good_per_frame.index(max(prop_good_per_frame)))
		sys.stdout.write(prog+"\r")
		sys.stdout.flush()
		alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
		if goodlocus == True:
			goodlocus = False
			if opt == "-t":
				fnew = f.split("/")
				fn = fnew[len(fnew)-1]
				fn2 = "./translated/"+".".join(fn.split(".")[:-1])+".fas"
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
					print (transseq, file=debug)
					print (t, file=debug)
					transseq = str(transseq)
					#if remainder > 0:
					for t in range(len(transseq)-1,-1,-1):
						if transseq[t] != "X":
							transseq = transseq[:t]+"X"+transseq[t+1:]
							break
					print (">"+seq.id, "\n", transseq, file=outfile)#unambig_translate
				outfile.close()
			elif opt == "-orf":
				print (f, "\t", frame, file=outfile)
			print ("Final frame:", frame, file=debug)
			print ("---------------------------------------------------------", file=debug)
		else:
			if opt == "-orf":
				print (f, "\t", frame, file=outfile)
			print ("PREDICTED final frame:", frame, file=debug)


			for frame_num in range(0,6):
				print ("PREDICTED alignment for frame", frame_num, file=debug)
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
					print (transseq, file=debug)
					#print >> debug, t
			#alignment.seek(0)
			print ("---------------------------------------------------------", file=debug)
elif opt == "-u":
	if not os.path.exists ("./nt_translator"):
		os.makedirs("./nt_translator")
	for f in files:
		prog = "working on file "+str(f)
		sys.stdout.write(prog+"\r")
		sys.stdout.flush()
		seqs = SeqIO.parse(f, "fasta")
		fnew = f.split("/")
		fn = fnew[len(fnew)-1]
		fn2 = "./translated/"+".".join(fn.split(".")[:-1])+".fas"
		fn3 = "./nt_translator/"+".".join(fn.split(".")[:-1])+".fas"
		outfile = open(fn2, "w")
		outfile2 = open(fn3, "w")
		for seq in seqs:
			transseq = str(seq.seq.translate()).upper()#.replace("*","X")
			stopindex = [pos for pos, char in enumerate(transseq) if char == "*"]
			if len(stopindex) > 0:
				filtered_nt = ""
				filtered_aa = ""
				for basenum in range(0,len(seq.seq),3):
					if int(basenum/3) not in stopindex:
						filtered_aa += transseq[int(basenum/3)]
						filtered_nt += seq.seq[basenum:basenum+3]
				if len(filtered_nt) > 0 and len(filtered_aa) > 0:
					print (">"+seq.id.rstrip(), file=outfile)
					print (filtered_aa, file=outfile)
					print (">"+seq.id.rstrip(), file=outfile2)
					print (filtered_nt, file=outfile2)
			else:
				print (">"+seq.id.rstrip(), file=outfile)
				print (transseq, file=outfile)
				print (">"+seq.id.rstrip(), file=outfile2)
				print (seq.seq, file=outfile2)
		outfile.close()
		outfile2.close()
if opt == "-orf":
	outfile.close()
if opt == "-t" or opt == "-orf":
	print ("")
	print ("good loci:", count)
	for x, y in badloci.items():
		print (x, "predicted frame:", y)
	print ("bad loci:", len(badloci))
else:
	print ("")
print ("done")