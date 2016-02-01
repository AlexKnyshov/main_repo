from Bio import AlignIO
import sys
import glob
import os
if len(sys.argv) == 7:
	inputfolder = sys.argv[1]
	inputformat = sys.argv[2]
	inputext = sys.argv[3]
	outputfolder = sys.argv[4]
	outputformat = sys.argv[5]
	outputext = sys.argv[6]
elif len(sys.argv) == 5:
	infilename = sys.argv[1]
	inputformat = sys.argv[2]
	outputformat = sys.argv[3]
	outfileext = sys.argv[4]
else:
	print "FORMAT: python fconv.py [inputfolder] [inputformat] [inputext] [outputfolder] [outputformat] [outputext]"
	print "EXAMPLE: python fconv.py ./fasta fasta fas ./phylip phylip-relaxed phylip"
	print "FORMAT: python fconv.py [inputfile] [inputformat] [outputformat] [outputext]"
	print "EXAMPLE: python fconv.py ./test.fas fasta phylip-relaxed .phy"
	sys.exit()

if len(sys.argv == 7):
	files = glob.glob(inputfolder+"/*."+inputext)
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
		alignments = AlignIO.parse(input_handle, inputformat)
		AlignIO.write(alignments, output_handle, outputformat)
		output_handle.close()
		input_handle.close()
		count += 1
	print "Converted %i records" % count
elif len(sys.argv == 5):
	fnew = infilename+outfileext
	input_handle = open(infilename, "rU")
	output_handle = open(fnew, "w")
	print "converting", infilename, "to", fnew
	alignments = AlignIO.parse(input_handle, inputformat)
	AlignIO.write(alignments, output_handle, outputformat)
	output_handle.close()
	input_handle.close()	
print "Done"