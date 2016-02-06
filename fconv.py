from Bio import AlignIO
from Bio import SeqIO
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
	outfileext = sys.argv[5]
else:
	print "fconv.py script for converting aligned or unaligned sequence files"
	print "-----------folder input------------"
	print "FORMAT: python fconv.py [option: -a, -s] [inputfolder] [inputformat] [inputext] [outputfolder] [outputformat] [outputext]"
	print "EXAMPLE: python fconv.py -a ./fasta fasta .fas ./phylip phylip-relaxed .phy"
	print "------------file input-------------"
	print "FORMAT: python fconv.py [option: -a, -s] [inputfile] [inputformat] [outputformat] [outputext]"
	print "EXAMPLE: python fconv.py -a ./test.fas fasta phylip-relaxed .phy"
	print "--------some format options--------"
	print "fasta - fasta format"
	print "phylip - basic phylip with truncated names"
	print "phylip-relaxed - extended phylip (only in -a mode)"
	sys.exit()

if len(sys.argv) == 8:
	files = glob.glob(inputfolder+"/*"+inputext)
	extlen = len(inputext)
	if not os.path.exists ("./"+outputfolder):
		os.makedirs("./"+outputfolder)
	count = 0
	for f in files:
		fnew = f.split("/")
		fn = fnew[len(fnew)-1]
		print fn[:(len(fn)-extlen)]
		print outputfolder+"/"+fn[:(len(fn)-extlen)]+outputext
		input_handle = open(f, "rU")
		output_handle = open(outputfolder+"/"+fn[:(len(fn)-extlen)]+outputext, "w")
		if option == "-a":
			alignments = AlignIO.parse(input_handle, inputformat)
			AlignIO.write(alignments, output_handle, outputformat)
		elif option == "-s":
			sequences = SeqIO.parse(input_handle, inputformat)
			SeqIO.write(sequences, output_handle, outputformat)
		output_handle.close()
		input_handle.close()
		count += 1
	print "Converted %i records" % count
elif len(sys.argv) == 6:
	fnew = infilename+outfileext
	input_handle = open(infilename, "rU")
	output_handle = open(fnew, "w")
	print "converting", infilename, "to", fnew
	if option == "-a":
		alignments = AlignIO.parse(input_handle, inputformat)
		AlignIO.write(alignments, output_handle, outputformat)
	elif option == "-s":
		sequences = SeqIO.parse(input_handle, inputformat)
		SeqIO.write(sequences, output_handle, outputformat)
	output_handle.close()
	input_handle.close()	
print "Done"