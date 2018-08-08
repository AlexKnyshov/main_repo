from Bio import SeqIO
import csv
import sys
import glob
import os
import shutil

if len(sys.argv) == 4:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*.fas")
	if sys.argv[2] == "-a" or sys.argv[2] == "-ar" or sys.argv[2] == "-e":
		exclusion_file = sys.argv[3]
	elif sys.argv[2] == "-l" or sys.argv[2] == "-ll":
		threshold = float(sys.argv[3])
	else:
		print "incorrect command line parameters"
		sys.exit()
else:
	print "FORMAT: python removeTaxa.py [folder with fasta] [option: -a (leave specified taxa), -ar (leave and rename specified taxa), -e (exclude specified taxa), -l (exclude short seq taxa), -ll (exclude loci with few taxa)] [taxalist (or csv) or lenght percent threshold]"
	print "EXAMPLE: python removeTaxa.py ./fasta -l 0.75"
	print "EXAMPLE: python removeTaxa.py ./fasta -a list.lst"
	sys.exit()

if sys.argv[2] == "-a" or sys.argv[2] == "-e":
	print "reading taxalist..."
	exclusion_list = []
	exfile = open(exclusion_file, "r")
	for line in exfile:
		l = line.strip()
		exclusion_list.append(l)
	exfile.close()
	print "read", len(exclusion_list), "records"

if sys.argv[2] == "-ar":
	print "reading csv taxalist..."
	exclusion_list = {}
	exfile = open(exclusion_file, "r")
	reader = csv.reader(exfile)
	for row in reader:
		#l = line.strip()
		exclusion_list[row[0]] = row[1]
	exfile.close()
	print "read", len(exclusion_list), "records"

print "creating an output folder..."
if not os.path.exists ("./reduced/"):
    os.makedirs("./reduced") #creating folder if necessary
else:
    shutil.rmtree("./reduced/") #removing old files
    os.makedirs("./reduced")

print "parsing the files..."
for f in files:
	fn = f.split("/")[-1]
	prog = "working on file "+fn
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
 	if sys.argv[2] == "-ll":
 		#print len([SeqIO.parse(f, "fasta")])
 		if len(list(SeqIO.parse(f, "fasta"))) >= threshold:
 			shutil.copy2(inputfolder+fn, "./reduced")
 	else:
		outputfile=open("./reduced/"+fn, "w")
		count = 0
		for seq in SeqIO.parse(f, "fasta"):
			#print seq.id, len(str(seq.seq).replace("-", "").replace("N", "")), len(seq.seq), float(len(str(seq.seq).replace("-", "").replace("N", "")))/len(seq.seq)
			if sys.argv[2] == "-e" and seq.id not in exclusion_list:
				print >> outputfile, ">"+seq.id, "\n", seq.seq
				count += 1
			elif sys.argv[2] == "-a" and seq.id in exclusion_list:
				print >> outputfile, ">"+seq.id, "\n", seq.seq
				count += 1
			elif sys.argv[2] == "-ar" and seq.id in exclusion_list:
				print >> outputfile, ">"+exclusion_list[seq.id], "\n", seq.seq
				count += 1
			elif sys.argv[2] == "-l" and float(len(str(seq.seq).replace("-", "").upper().replace("N", "")))/len(seq.seq)>threshold:
				print >> outputfile, ">"+seq.id, "\n", seq.seq
				count += 1
		outputfile.close()
		if count == 0:
			os.remove("./reduced/"+fn)
print "\ndone"
