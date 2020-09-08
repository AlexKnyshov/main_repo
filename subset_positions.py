import os
import sys
import glob
import shutil
from Bio import AlignIO

#filepath input
if len(sys.argv) == 2:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*")
else:
	print "FORMAT: python subset_positions.py [folder with files]"
	print "EXAMPLE: python subset_positions.py ./folder"
	sys.exit()
if len(files) == 0:
	print "no files in the directory"


print "creating an output folder..."
if not os.path.exists ("./subset/"):
	os.makedirs("./subset") #creating folder if necessary
else:
	shutil.rmtree("./subset/") #removing old files
	os.makedirs("./subset")


for f in files:
	infile = open(f, "r")
	inalignment = AlignIO.read(infile, "fasta")
	outalignment = inalignment[:, :2]
	for pos in range(3, inalignment.get_alignment_length(),3):
		outalignment += inalignment[:, pos:pos+2]
	outpath = "./subset/"+f.split("/")[-1]
	AlignIO.write(outalignment, outpath, "fasta")
print "done"