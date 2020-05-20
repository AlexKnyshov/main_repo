from Bio import SeqIO
import csv
import sys
import glob
import os
import shutil

if len(sys.argv) == 3:
	csvfilename = sys.argv[1]
	seqfolder = sys.argv[2]
else:
	print "FORMAT: python UPhO_parser.py [csv file] [folder with alignments]"
	print "EXAMPLE: python UPhO_parser.py UPhO_nr_orthogroups.csv ./fasta/"
	sys.exit()

print "creating an output folder..."
if not os.path.exists ("./reduced/"):
    os.makedirs("./reduced") #creating folder if necessary
else:
    shutil.rmtree("./reduced/") #removing old files
    os.makedirs("./reduced")

maindict = {}
csvfile = open(csvfilename, "r")
reader = csv.reader(csvfile)
for row in reader:
	locus = "_".join(row[0].split("/")[1].split("_")[:-1])
	maindict[locus] = []
	for r in row[1:]:
		if r != '':
			maindict[locus].append(r)
csvfile.close()
print maindict
sys.exit()
for key, value in maindict.items():
	if len(maindict) == 1:
		orig = open(seqfolder+"/"+f.split("bipartitions.")[1].split("_UPhO")[0])
		outputfile=open("./reduced/"+f.split("bipartitions.")[1].split("_UPhO")[0], "w")
		for seq in SeqIO.parse(orig, "fasta"):
			if seq.id in value:
				print >> outputfile, ">"+seq.id.split("|")[0], "\n", seq.seq
	else:
		print "skip for now"
		#len(value)
#print "read", len(exclusion_list), "records"
outputfile.close()
orig.close()
