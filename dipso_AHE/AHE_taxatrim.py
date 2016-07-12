from Bio import SeqIO
import sys
import glob
import os
import shutil

if len(sys.argv) == 3:
    f = sys.argv[1]
    files = glob.glob(f+"/*")
    listfile = sys.argv[2]
else:
    print "FORMAT: python AHE_taxatrim.py [folder with fasta] [trimlist]"
    print "EXAMPLE: python AHE_taxatrim.py ./fasta list.lst"
    sys.exit()

translist = []
if len(sys.argv) == 6:
    print "reading exclusion list..."
    lfile = open(listfile, "r")
    for line in lfile:
        l = line.strip()
        translist.append(l)
    lfile.close()
    print "read", len(translist), "records"

locilist = set()
print "input"
for infile in files:
	input_handle = open(infile, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	print "read", infile
	#outhandle = open(infile+".tx_tm.fas", "w")
	for seq in alignments:
		#if seq.id in trimlist:
		if seq.id in translist:
			locilist.add(infile)
#print locilist
print len(locilist)

#copy files
if not os.path.exists ("./subset/"):
    os.makedirs("./subset") #creating folder if necessary
else:
    shutil.rmtree("./subset/") #removing old files
    os.makedirs("./subset")

print "copying files:"
for x in locilist:
    locusfname = x.split("/")[-1]
    #print locusfname
    if not os.path.exists ("./subset/"+locusfname):
        prog = "copying "+str(locusfname)+"..."
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        shutil.copy2(f+locusfname, "./subset")
print "done"
