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
	print ("FORMAT: python customtrim.py [folder with fasta] [trimming option: -a, -1, -%, -refine, -d, -d%, -tiger] ([trimlist])")
	print ("EXAMPLE: python customtrim.py ./fasta -%")
	print ("EXAMPLE: python customtrim.py ./fasta -1 trimlist.txt")
	print ("EXAMPLE: python customtrim.py ./fasta -tiger .tig")
	sys.exit()
if len(files) == 0:
	print ("no fasta files in the directory")

print ("the option", trimopt, "set up")
if trimopt[:2] == "-%":
	print ("the cutoff is set to", trimopt[2:])
elif trimopt[:3] == "-d%":
	print ("the cutoff is set to", trimopt[3:])
#starting to process files
loci = {}
if trimopt == "-a" or trimopt == "-1" or trimopt[:3] == "-d%":
	print ("reading taxalist...")
	taxalist = []
	exfile = open(exclusion_file, "r")
	for line in exfile:
		l = line.strip()
		taxalist.append(l)
	exfile.close()
	print ("read", len(taxalist), "records")
warninglist = []
progbarc = 0
print ("creating log file and output folder...")
outf = open("customtrim.out", "w")
print ("command line parameters:",' '.join(sys.argv), file=outf)

if not os.path.exists ("./trimmed"):
	os.makedirs("./trimmed")

print ("parsing files:")
for f in files:
	reftaxa = False
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	fn2 = "./trimmed/"+fnew[-1]#fn.split(".")[0]+".fas"
	print ("input:", f, file=outf)
	print ("output:", fn2, file=outf)
	infile = open(f, "r")
	seqs = {}
	inputalignment = AlignIO.read(infile, "fasta")
	for seq in inputalignment:
	    seqs[seq.id] = str(seq.seq).upper()
	    if trimopt == "-a" or trimopt == "-1" or trimopt[:3] == "-d%":
		    if seq.id in taxalist:
		    	reftaxa = True
	print ("reftaxa:", reftaxa, file=outf)
	#pairwise
	names =[]
	lengths = []
	if trimopt[:2] == "-%" or trimopt == "-refine" or trimopt == "-d" or reftaxa == False: #trim counting all taxa
		for key in seqs.keys():
			names.append(key)
			lengths.append(len(seqs[key]))
	elif reftaxa: #trim to AHE taxa
		for key in seqs.keys():
			if key in taxalist:
				names.append(key)
				lengths.append(len(seqs[key]))
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
				break
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
				datablock += 1 #let's start counting...
				dataflag = True #full block
				if datablock == 1:
					potential_startpos = startpos #only remember this when we entered the block
			else: #gappy position
				if dataflag == True: #already in block
					if datablock > 20: #check how large was the block
						break
					else:
						datablock = 0
						dataflag = False
			startpos += 1
		print ("datablock (L to R): ", datablock, file=outf)
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
		positions = []
		for basenum in range(len(seqs[names[0]])):
			if basenum+1 not in badpos:
				positions.append(basenum)
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
				break
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
				break
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
				datablock += 1 #let's start counting...
				dataflag = True #full block
				if datablock == 1:
					potential_endpos = endpos #only remember this when we entered the block
			else: #gappy position
				if dataflag == True: #already in block
					if datablock > 20: #check how large was the block
						break
					else:
						datablock = 0
						dataflag = False
			endpos -= 1
		print ("datablock (R to L): ", datablock, file=outf)
		endpos = potential_endpos
	infile.close()
	if trimopt == "-tiger" or trimopt == "-d" or trimopt[:3] == "-d%" and reftaxa == True:
		print (len(positions),"good positions:", positions, file=outf)
	else:
		print ("startpos: ", startpos, file=outf)
		print ("endpos: ", endpos, file=outf)

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
			for aliseq in align:
				print (">"+aliseq.id, file=outfile)
				print (aliseq.seq, file=outfile)
			outfile.close()
			for aliseq in alignmasked:
				print (">"+aliseq.id, file=outmask)
				print (aliseq.seq, file=outmask)
			outmask.close()
		elif (trimopt == "-tiger" or trimopt == "-d" or trimopt[:3] == "-d%" and reftaxa == True) and len(positions) == 0:
			print (len(positions), "good positions, skipping the locus...", file=outf)
			warninglist.append(f)
		else:
			outfile = open(fn2, "w")
			outmask = open(fn2+"_masked", "w")
			for tempseq in inputalignment:
				print (">"+tempseq.id, file=outfile)
				print (tempseq.seq[startpos:endpos], file=outfile)
				print (">"+tempseq.id, file=outmask)
				print (("$"*startpos)+tempseq.seq[startpos:endpos]+("$"*(len(tempseq.seq)-endpos-1)), file=outmask)
			outfile.close()
			outmask.close()
	else:
		warninglist.append(f)
	

	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	if trimopt == "-d" or trimopt[:3] == "-d%" or trimopt == "-tiger":
		prog = str(progbar)+"% working on file "+str(f)+": "+str(len(positions))+" good positions"
	else:
		prog = str(progbar)+"% working on file "+str(f)+": starts "+str(startpos)+", ends "+str(endpos)
	sys.stdout.write(prog+"\r")
	sys.stdout.flush()
print ("\nwarning list:", len(warninglist))
if len(warninglist)>0:
	for x in warninglist:
		print (x)
outf.close()
print ("done")