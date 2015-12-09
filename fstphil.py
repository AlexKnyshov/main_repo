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
else:
	print "FORMAT: python fstphil.py [inputfolder] [inputformat] [inputext] [outputfolder] [outputformat] [outputext]"
	print "EXAMPLE: python fstphil.py ./fasta fasta fas ./phylip phylip-relaxed phylip"
	sys.exit()

files = glob.glob(inputfolder+"/*."+inputext)
extlen = len(inputext)
if not os.path.exists ("./"+outputfolder):
	os.makedirs("./"+outputfolder)
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
	#count = SeqIO.convert(f, inputformat, outputfolder+"/"+fn[:(len(fn)-extlen)]+"phylip", outputformat)
#print "Converted %i records" % count
print "Done"