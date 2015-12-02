from Bio import SeqIO
import sys
inputfolder = sys.argv[1]
inputformat = sys.argv[2]
outputfolder = sys.argv[3]
outputformat = sys.argv[4]
count = SeqIO.convert(inputfolder, inputformat, outputfolder, outputformat)
print "Converted %i records" % count
print "Done"