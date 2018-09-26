from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
import sys
import glob
import os
if len(sys.argv) == 8:
	option = sys.argv[1]
	inputfolder = sys.argv[2]
	inputformat = sys.argv[3]
	inputext = sys.argv[4]
	outputfolder = sys.argv[5]
	outputformat = sys.argv[6]
	outputext = sys.argv[7]
elif len(sys.argv) == 6:
	option = sys.argv[1]
	infilename = sys.argv[2]
	inputformat = sys.argv[3]
	outputformat = sys.argv[4]
	outputext = sys.argv[5]
else:
	print "fconv.py script for converting aligned or unaligned sequence files"
	print "-----------folder input------------"
	print "FORMAT: python fconv.py [option: -a, -s, -print] [inputfolder] [inputformat] [inputext] [outputfolder] [outputformat] [outputext]"
	print "EXAMPLE: python fconv.py -a ./fasta fasta .fas ./phylip phylip-relaxed .phy"
	print "------------file input-------------"
	print "FORMAT: python fconv.py [option: -a, -s, -print] [inputfile] [inputformat] [outputformat] [outputext]"
	print "EXAMPLE: python fconv.py -a ./test.fas fasta phylip-relaxed .phy"
	print "--------general manual--------"
	print "for options -a (AlignIO) and -s (SeqIO) see format options in the corresponding biopython module manual. Some are listed below"
	print "fasta - fasta format"
	print "phylip - basic phylip with truncated names"
	print "phylip-relaxed - extended interleaved phylip (only in -a mode)"
	print "option -print in conjunction with output format set to fasta writes fasta file directly"
	print "option -print in conjunction with output format set to phylip-relaxed writes non-interleaved phylip-relaxed file directly"
	sys.exit()

def check_alphabet(filename, fformat):#code for this function is modified from https://stackoverflow.com/questions/41588725/estimate-alphabet-in-biopython-from-fasta-file/41596563#41596563
	alphabets = [IUPAC.ambiguous_dna, IUPAC.protein]#[ambiguous_dna, extended_protein]#, unambiguous_dna, extended_dna, ambiguous_rna, unambiguous_rna]
	first_record = list(SeqIO.parse(filename, fformat))[0]
	detected = ""
	#check DNA first:
	leftover = set(str(first_record.seq).upper()) - set(alphabets[0].letters) - set(["-", "?"])
	if len(leftover) == 0:
		detected = "DNA"
	else:
		leftover = set(str(first_record.seq).upper()) - set(alphabets[1].letters) - set(["-", "?","X"])
		if len(leftover) == 0:
			detected = "Prot"
		else:
			print filename, "error: unknown alphabet, problematic symbols:", leftover
			sys.exit()
	return detected

if len(sys.argv) == 8:
	files = glob.glob(inputfolder+"/*"+inputext)
	extlen = len(inputext)
	if not os.path.exists ("./"+outputfolder):
		os.makedirs("./"+outputfolder)
else:
	files = [infilename]
	extlen = 0
	outputfolder = "."

count = 0
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	print fn[:(len(fn)-extlen)]
	print outputfolder+"/"+fn[:(len(fn)-extlen)]+outputext
	input_handle = open(f, "rU")
	output_handle = open(outputfolder+"/"+fn[:(len(fn)-extlen)]+outputext, "w")
	if option == "-a":
		if outputformat == "nexus" or inputformat == "nexus":
			alph = check_alphabet(input_handle, inputformat)
			input_handle.seek(0)
			print alph
			if alph == "DNA":
				alignments = AlignIO.parse(input_handle, inputformat, alphabet = Gapped(IUPAC.ambiguous_dna))
			elif alph == "Prot":
				alignments = AlignIO.parse(input_handle, inputformat, alphabet = Gapped(IUPAC.protein, '-'))
		else:
			alignments = AlignIO.parse(input_handle, inputformat)
		#print list(alignments)
		AlignIO.write(alignments, output_handle, outputformat)
	elif option == "-s":
		if outputformat == "nexus" or inputformat == "nexus":
			alph = check_alphabet(input_handle, inputformat)
			input_handle.seek(0)
			print alph
			if alph == "DNA":
				sequences = SeqIO.parse(input_handle, inputformat, alphabet = Gapped(IUPAC.ambiguous_dna))
			elif alph == "Prot":
				sequences = SeqIO.parse(input_handle, inputformat, alphabet = Gapped(IUPAC.protein, '-'))
		else:
			sequences = SeqIO.parse(input_handle, inputformat)
		SeqIO.write(sequences, output_handle, outputformat)
	elif option == "-print" and outputformat == "fasta":
		sequences = SeqIO.parse(input_handle, inputformat)
		for seq in sequences:
			print >> output_handle, ">"+seq.id, "\n", seq.seq
	elif option == "-print" and outputformat == "phylip-relaxed":
		alignments = AlignIO.read(input_handle, inputformat)
		print >> output_handle, str(len(alignments))+" "+str(alignments.get_alignment_length())
		for seq in alignments:
			print >> output_handle, str(seq.id)+" "+str(seq.seq)
	output_handle.close()
	input_handle.close()
	count += 1
print "Converted %i records" % count
print "Done"