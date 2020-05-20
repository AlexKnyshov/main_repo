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
	maindict[row[0]] = []
	for r in row[1:]:
		if r != '':
			maindict[row[0]].append(r)
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
# print "parsing the files..."
# for f in files:
# 	fn = f.split("/")[-1]
# 	prog = "working on file "+fn
#  	sys.stdout.write(prog+"\r")
#  	sys.stdout.flush()
#  	if sys.argv[2] == "-ll":
#  		#print len([SeqIO.parse(f, "fasta")])
#  		if len(list(SeqIO.parse(f, "fasta"))) >= threshold:
#  			shutil.copy2(inputfolder+fn, "./reduced")
#  	else:
# 		outputfile=open("./reduced/"+fn, "w")
# 		count = 0
# 		for seq in SeqIO.parse(f, "fasta"):
# 			#print seq.id, len(str(seq.seq).replace("-", "").replace("N", "")), len(seq.seq), float(len(str(seq.seq).replace("-", "").replace("N", "")))/len(seq.seq)
# 			if sys.argv[2] == "-e" and seq.id not in exclusion_list:
# 				print >> outputfile, ">"+seq.id, "\n", seq.seq
# 				count += 1
# 			elif sys.argv[2] == "-a" and seq.id in exclusion_list:
# 				print >> outputfile, ">"+seq.id, "\n", seq.seq
# 				count += 1
# 			elif sys.argv[2] == "-ar" and seq.id in exclusion_list:
# 				print >> outputfile, ">"+exclusion_list[seq.id], "\n", seq.seq
# 				count += 1
# 			elif sys.argv[2] == "-l" and float(len(str(seq.seq).replace("-", "").upper().replace("N", "").replace("?", "")))/len(seq.seq)>threshold:
# 				print >> outputfile, ">"+seq.id, "\n", seq.seq
# 				count += 1
# 		outputfile.close()
# 		if count == 0:
# 			os.remove("./reduced/"+fn)
# print "\ndone"
