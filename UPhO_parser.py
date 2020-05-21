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
	taxa = []
	for r in row[1:]:
		if r != '':
			taxa.append(r)
	if locus in maindict:
		if len(taxa) > len(maindict[locus]):
			maindict[locus] = taxa
	else:
		maindict[locus] = taxa
	
csvfile.close()

for key, value in maindict.items():
	orig = open(seqfolder+"/"+key+".fas")
	print "processing", key
	outputfile=open("./reduced/"+key+".fas", "w")
	for seq in SeqIO.parse(orig, "fasta"):
		if seq.id in value:
			print >> outputfile, ">"+seq.id.split("|")[0], "\n", seq.seq
	outputfile.close()
	orig.close()