from Bio import SeqIO
import sys
import glob
import os
inputfolder = sys.argv[1]
inputformat = sys.argv[2]
inputext = sys.argv[3]
outputfolder = sys.argv[4]
outputformat = sys.argv[5]
files = glob.glob(inputfolder+"/*."+inputext)
extlen = len(inputext)
if not os.path.exists ("./"+outputfolder):
	os.makedirs("./"+outputfolder)
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	print fn[:(len(fn)-extlen)]
	print outputfolder+"/"+fn[:(len(fn)-extlen)]+"phylip"
	count = SeqIO.convert(f, inputformat, outputfolder+"/"+fn[:(len(fn)-extlen)]+"phylip", outputformat)
#print "Converted %i records" % count
print "Done"