import os
import sys
import glob
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
#filepath input
if len(sys.argv) >= 3:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*.fas")
	trimopt = sys.argv[2]
	if trimopt == "-a" or trimopt == "-1" or trimopt[:3] == "-d%":
		exclusion_file = sys.argv[3]
	elif trimopt == "-tiger":
		tiger_file = sys.argv[3]
	elif trimopt == "-d":
		ff = sys.argv[3]
		if ff == "dna":
			alph = Gapped(IUPAC.ambiguous_dna)
		elif ff == "prot":
			alph = Gapped(IUPAC.protein, '-')
else:
	print "FORMAT: python customtrim.py [folder with fasta] [trimming option: -a, -1, -%, -refine, -d, -d%, -tiger] ([trimlist])"
	print "EXAMPLE: python customtrim.py ./fasta -%"
	print "EXAMPLE: python customtrim.py ./fasta -1 trimlist.txt"
	print "EXAMPLE: python customtrim.py ./fasta -tiger .tig"
	sys.exit()
if len(files) == 0:
	print "no fasta files in the directory"

print "the option", trimopt, "set up"
if trimopt[:2] == "-%":
	print "the cutoff is set to", trimopt[2:]
elif trimopt[:3] == "-d%":
	print "the cutoff is set to", trimopt[3:]
#starting to process files
print "initializing..."
loci = {}
#taxalist=["I13432_ED_4993_Hemiptera_Dipsocoridae_Cryptostemma_sp_seq1", "I13433_ED_2045_Hemiptera_Ceratocombidae_Kvamula_sp_seq1", "I13434_ED_2660_Hemiptera_Schizopteridae_Williamsocoris_sp_seq1", "I13435_ED_4258_Hemiptera_Schizopteridae_Nannocoris_sp_seq1", "I13436_ED_1692_Hemiptera_Schizopteridae_Kokeshia_sp_seq1", "I13437_ED_4257_Hemiptera_Schizopteridae_Chinannus_sp_seq1", "I13438_ED_2192_Hemiptera_Schizopteridae_Dundonannus_sp_seq1", "I13439_ED_6303_Hemiptera_Schizopteridae_Schizoptera_sp_seq1", "I13440_P14_RCW_1261_Hemiptera_Reduviinae_Opisthacidius_sp_seq1", "I13442_UCRC_ENT_00092725_Hemiptera_Tribelocephalinae_Afrodecius_sp_seq1", "I13443_RCW_4586_Hemiptera_Vesciinae_Vescia_sp_seq1", "I13444_RCW_4525_Hemiptera_Reduviinae_Rulandus_phaedrus_seq1", "I13445_RCW_4101_Hemiptera_Phymatinae_Phymata_pacifica_seq1"]
#taxalist=["Chinannus_monteverdensis_4257", "Schizoptera_sp_6303", "cf_Kvamula_or_Seychellesanus_sp_Madagascar_2043", "Williamsocoris_sp_Trinidad_2660", "Nannocoris_sp_4258", "Dundonannus_sp_2190", "Cryptostemma_sp_Peru_249", "Kokeshia_sp_Thailand_1409"]
if trimopt == "-a" or trimopt == "-1" or trimopt[:3] == "-d%":
	print "reading taxalist..."
	taxalist = []
	exfile = open(exclusion_file, "r")
	for line in exfile:
		l = line.strip()
		taxalist.append(l)
	exfile.close()
	print "read", len(taxalist), "records"
warninglist = []
progbarc = 0
print "creating log file and output folder..."
outf = open("customtrim.out", "w")
print >> outf, "command line parameters:",' '.join(sys.argv)

if not os.path.exists ("./trimmed"):
	os.makedirs("./trimmed")

print "parsing files:"
for f in files:
	reftaxa = False
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	fn2 = "./trimmed/"+fnew[-1]#fn.split(".")[0]+".fas"
	print >> outf, "input:", f
	print >> outf, "output:", fn2
	infile = open(f, "r")
	seqs = {}
	inputalignment = AlignIO.read(infile, "fasta")
	for seq in inputalignment:
	    seqs[seq.id] = str(seq.seq).upper()
	    if trimopt == "-a" or trimopt == "-1" or trimopt[:3] == "-d%":
		    if seq.id in taxalist:
		    	reftaxa = True
	print >> outf, "reftaxa:", reftaxa
	#pairwise
	names =[]
	lengths = []
	if trimopt[:2] == "-%" or trimopt == "-refine" or trimopt == "-d" or reftaxa == False: #trim counting all taxa
		for key in seqs.keys():
			names.append(key)
			lengths.append(len(seqs[key]))
			#print key, len(seqs[key])	
	elif reftaxa: #trim to AHE taxa
		for key in seqs.keys():
			if key in taxalist:
				names.append(key)
				lengths.append(len(seqs[key]))
				#print key, len(seqs[key])
	#else skip	
	distlist = []
	#forward trim
	startpos = 0
	endpos = max(lengths)
	if trimopt == "-a" and reftaxa == True:
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
	elif trimopt == "-1" and reftaxa == True:
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
	elif trimopt[:2] == "-%":
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata < (len(names)*float(trimopt[2:])):
				break
			#print misdata, "misdata", startpos, "startpos"
			startpos += 1
	elif trimopt == "-refine":
		datablock = 0
		dataflag = False
		potential_startpos = 0
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1 #mis data in a position
			if misdata <= 0: #full position
				#print "done"
				datablock += 1 #let's start counting...
				dataflag = True #full block
				if datablock == 1:
					potential_startpos = startpos #only remember this when we entered the block
			else: #gappy position
				if dataflag == True: #already in block
					if datablock > 20: #check how large was the block
						#print >> outf, "datablock (L to R): ", datablock
						#startpos = potential_startpos
						break
					else:
						datablock = 0
						dataflag = False
			#print misdata, "misdata", startpos, "startpos"
			startpos += 1
		print >> outf, "datablock (L to R): ", datablock
		startpos = potential_startpos
	elif trimopt == "-d":
		positions = []
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "?":
					misdata += 1
				if base == "N" and ff == "dna":
					misdata += 1
				if base == "X" and ff == "prot":
					misdata += 1
				if misdata > (len(names)*0.4):
					break
			if misdata <= (len(names)*0.4):
				positions.append(basenum)
	elif trimopt[:3] == "-d%" and reftaxa == True:
		positions = []
		for basenum in range(len(seqs[names[0]])):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
				if misdata > (len(names)*float(trimopt[3:])):
					break
			if misdata <= (len(names)*float(trimopt[3:])):
				positions.append(basenum)
	elif trimopt == "-tiger":
		with open(f+tiger_file) as tf:
			badpos = [int(s) for s in tf.readline().strip().split("=")[1].strip(";").split()]
			#badpos = tf.readline().strip().split("=")[1].strip(";").split()
		positions = []
		for basenum in range(len(seqs[names[0]])):
			if basenum+1 not in badpos:
				positions.append(basenum)
	#print startpos, "startpos"
	#reverse trim
	if trimopt == "-a" and reftaxa == True:
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
	elif trimopt == "-1" and reftaxa == True:
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
	elif trimopt[:2] == "-%":
		for basenum in range(len(seqs[names[0]])-1, -1, -1):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata < (len(names)*float(trimopt[2:])):
				break
			#print misdata, "misdata", startpos, "startpos"
			endpos -= 1
	elif trimopt == "-refine":
		datablock = 0
		dataflag = False
		potential_endpos = endpos
		for basenum in range(len(seqs[names[0]])-1, -1, -1):
			baselist = []
			misdata = 0
			for name in names:
				baselist.append(seqs[name][basenum])
			for base in baselist:
				if base == "-" or base == "N" or base == "?":
					misdata += 1
			if misdata <= 0: #full position
				#print "done"
				datablock += 1 #let's start counting...
				dataflag = True #full block
				if datablock == 1:
					potential_endpos = endpos #only remember this when we entered the block
			else: #gappy position
				if dataflag == True: #already in block
					if datablock > 20: #check how large was the block
						#print >> outf, "datablock (R to L): ", datablock
						#endpos = potential_endpos
						break
					else:
						datablock = 0
						dataflag = False
			#print misdata, "misdata", endpos, "endpos"
			endpos -= 1
		print >> outf, "datablock (R to L): ", datablock
		endpos = potential_endpos
	#print endpos, "endpos"
	infile.close()
	if trimopt == "-tiger" or trimopt == "-d" or trimopt[:3] == "-d%" and reftaxa == True:
		print >> outf, len(positions),"good positions:", positions
	else:
		print >> outf, "startpos: ", startpos
		print >> outf, "endpos: ", endpos

	#writing to files
	if endpos-startpos > 0 or trimopt == "-d" or trimopt[:3] == "-d%" and reftaxa == True:
		if (trimopt == "-tiger" or trimopt == "-d" or trimopt[:3] == "-d%" and reftaxa == True) and len(positions) > 0:
			###########################################
			###change code to output masked file too###
			outfile = open(fn2, "w")
			outmask = open(fn2+"_masked", "w")
			align = MultipleSeqAlignment([], alph) # new ali
			alignmasked = MultipleSeqAlignment([], alph) # new ali
			inpalilen = inputalignment.get_alignment_length()
			for tempseq in inputalignment:
				seqrec = SeqRecord(Seq("", alph), id=tempseq.id) #add taxa
				maskedseqrec = SeqRecord(Seq("", alph), id=tempseq.id) #add taxa
				for p in range(inpalilen):
					# good position
					if p in positions:
						seqrec += tempseq[p]
						maskedseqrec += tempseq[p]
					#bad position
					else:
						maskedseqrec += "$"
				if len(str(seqrec.seq).replace("-","")) / float(len(seqrec.seq)) >= 0.1:
					align.extend([seqrec])
					alignmasked.extend([maskedseqrec])
			#AlignIO.write(align, outfile, "fasta")
			for aliseq in align:
				print >> outfile, ">"+aliseq.id
				print >> outfile, aliseq.seq
			outfile.close()
			for aliseq in alignmasked:
				print >> outmask, ">"+aliseq.id
				print >> outmask, aliseq.seq
			outmask.close()
		elif (trimopt == "-tiger" or trimopt == "-d" or trimopt[:3] == "-d%" and reftaxa == True) and len(positions) == 0:
			print >> outf, len(positions), "good positions, skipping the locus..."
			warninglist.append(f)
		else:
			outfile = open(fn2, "w")
			outmask = open(fn2+"_masked", "w")
			for tempseq in inputalignment:
				print >> outfile, ">"+tempseq.id, "\n", tempseq.seq[startpos:endpos]
				print >> outfile, ">"+tempseq.id, "\n", ("$"*startpos)+tempseq.seq[startpos:endpos]+("$"*(len(tempseq.seq)-endpos-1))
			# for seq, s in seqs.items():
			# 	print >> outfile, ">"+seq, "\n", s[startpos:endpos]
			outfile.close()
			outmask.close()
	else:
		warninglist.append(f)
	

	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	#hashes = '#' * int(progbar * 0.2)
	#spaces = ' ' * (20 - len(hashes))
	#print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
	if trimopt == "-d" or trimopt[:3] == "-d%" or trimopt == "-tiger":
		prog = str(progbar)+"% working on file "+str(f)+": "+str(len(positions))+" good positions"
	else:
		prog = str(progbar)+"% working on file "+str(f)+": starts "+str(startpos)+", ends "+str(endpos)
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
print "\nwarning list:", len(warninglist)
if len(warninglist)>0:
	for x in warninglist:
		print x
outf.close()
print "done"