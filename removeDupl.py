from Bio import SeqIO
import sys
import glob
import os
import shutil

if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	files = glob.glob(inputfolder+"/*.fas")
	if sys.argv[2] == "-remove":
		remove = True
	elif sys.argv[2] == "-longest":
		remove = False
	else:
		print "incorrect command line parameters"
		sys.exit()
else:
	print "FORMAT: python removeDupl.py [folder with fasta] [option: -remove (remove taxa with duplicates), -longest (leave the longest seqeunce for each taxon)"
	print "EXAMPLE: python removeDupl.py ./fasta -remove"
	sys.exit()

print "creating an output folder..."
if not os.path.exists ("./rmDup/"):
    os.makedirs("./rmDup") #creating folder if necessary
else:
    shutil.rmtree("./rmDup/") #removing old files
    os.makedirs("./rmDup")

print "parsing the files..."
for f in files:
	fn = f.split("/")[-1]
	prog = "working on file "+fn
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
	outputfile=open("./rmDup/"+fn, "w")
	if remove: #remove option
		rmdict = {}
		for seq in SeqIO.parse(f, "fasta"):
			if seq.id in rmdict:
				rmdict[seq.id] += 1
			else:
				rmdict[seq.id] = 1
		#f.seek(0)
		for seq in SeqIO.parse(f, "fasta"):
			if rmdict[seq.id] == 1:
				print >> outputfile, ">"+seq.id, "\n", seq.seq
	else:
		rmdict = {}
		for seq in SeqIO.parse(f, "fasta"):
			if seq.id in rmdict:
				if len(str(rmdict[seq.id]).replace("-", "")) < len(str(seq.seq).replace("-", "")):
					rmdict[seq.id] = seq.seq
			else:
				rmdict[seq.id] = seq.seq
		for key, value in sorted(rmdict.items()):
			print >> outputfile, ">"+key, "\n", value
	outputfile.close()
print "\ndone"
