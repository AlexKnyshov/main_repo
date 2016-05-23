import os
import sys
import glob
from Bio import AlignIO
#filepath input
if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*.fas")
	trimopt = sys.argv[2]
else:
	print "FORMAT: python customtrim.py [folder with fasta] [trimming option: -a, -1, -%]"
	print "EXAMPLE: python customtrim.py ./fasta -1"
	sys.exit()
if len(files) == 0:
	print "no fasta files in the directory"

#starting to process files
print "initializing..."
loci = {}
taxalist=["I13432_ED_4993_Hemiptera_Dipsocoridae_Cryptostemma_sp_seq1", "I13433_ED_2045_Hemiptera_Ceratocombidae_Kvamula_sp_seq1", "I13434_ED_2660_Hemiptera_Schizopteridae_Williamsocoris_sp_seq1", "I13435_ED_4258_Hemiptera_Schizopteridae_Nannocoris_sp_seq1", "I13436_ED_1692_Hemiptera_Schizopteridae_Kokeshia_sp_seq1", "I13437_ED_4257_Hemiptera_Schizopteridae_Chinannus_sp_seq1", "I13438_ED_2192_Hemiptera_Schizopteridae_Dundonannus_sp_seq1", "I13439_ED_6303_Hemiptera_Schizopteridae_Schizoptera_sp_seq1", "I13440_P14_RCW_1261_Hemiptera_Reduviinae_Opisthacidius_sp_seq1", "I13442_UCRC_ENT_00092725_Hemiptera_Tribelocephalinae_Afrodecius_sp_seq1", "I13443_RCW_4586_Hemiptera_Vesciinae_Vescia_sp_seq1", "I13444_RCW_4525_Hemiptera_Reduviinae_Rulandus_phaedrus_seq1", "I13445_RCW_4101_Hemiptera_Phymatinae_Phymata_pacifica_seq1"]
warninglist = []
progbarc = 0
print "creating log file and output folder..."
outf = open("customtrim.out", "w")
if not os.path.exists ("./trimmed"):
	os.makedirs("./trimmed")
print "parsing files:"
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	fn2 = "./trimmed/"+fn.split(".")[0]+".fas"
	#print "input:", f
	#print "output:", fn2
	infile = open(f, "r")
	seqs = {}
	for seq in AlignIO.read(infile, "fasta"):
	    seqs[seq.id] = str(seq.seq).upper()
	#pairwise
	names =[]
	lengths = []
	#if trimopt == "-%":
	for key in seqs.keys():
		names.append(key)
		lengths.append(len(seqs[key])) #replace with a simplier solution - no need for this
		#print key, len(seqs[key])	
	# else:
	# 	for key in seqs.keys():
	# 		if key in taxalist:
	# 			names.append(key)
	# 			lengths.append(len(seqs[key])) #replace with a simplier solution - no need for this
	# 			#print key, len(seqs[key])	
	distlist = []
	#forward trim
	startpos = 0
	endpos = max(lengths)
	if trimopt == "-a":
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata <= 0:
				#print "done"
				break
			#print misdata, "misdata", startpos, "startpos"
			startpos += 1
	elif trimopt == "-1":
		finish = False
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
				else:
					finish = True
			if finish:
				break
			#print misdata, "misdata", startpos, "startpos"
			startpos += 1
	elif trimopt == "-%":
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata < (len(names)*0.1):
				break
			#print misdata, "misdata", startpos, "startpos"
			startpos += 1
	#print startpos, "startpos"
	#reverse trim
	if trimopt == "-a":
		for basenum in range(len(seqs[names[0]])-1, -1, -1):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata <= 0:
				#print "done"
				break
			#print misdata, "misdata", endpos, "endpos"
			endpos -= 1
	elif trimopt == "-1":
		finish = False
		for basenum in range(len(seqs[names[0]])-1, -1, -1):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
				else:
					finish = True
			if finish:
				#print "done"
				break
			#print misdata, "misdata", endpos, "endpos"
			endpos -= 1
	elif trimopt == "-%":
		for basenum in range(len(seqs[names[0]])-1, -1, -1):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata < (len(names)*0.1):
				break
			#print misdata, "misdata", startpos, "startpos"
			endpos -= 1
	#print endpos, "endpos"
	infile.close()
	print >> outf, "input:", f
	print >> outf, "output:", fn2
	print >> outf, startpos
	print >> outf, endpos
	if endpos-startpos > 0:
		outfile = open(fn2, "w")
		for seq, s in seqs.items():
			print >> outfile, ">"+seq, "\n", s[startpos:endpos]
		outfile.close()
	else:
		warninglist.append(f)
	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	#hashes = '#' * int(progbar * 0.2)
	#spaces = ' ' * (20 - len(hashes))
	#print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
	prog = str(progbar)+"% working on file "+str(f)+": starts "+str(startpos)+", ends "+str(endpos)
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
print "\nwarning list:", sum(warninglist)
outf.close()
print "done"