from Bio import SeqIO
import sys
import glob
import os
import shutil

if len(sys.argv) == 4:
    option = sys.argv[1]
    f = sys.argv[2]
    files = glob.glob(f+"/*")
    if option == "-t":
        listfile = sys.argv[3]
    elif option == "-l":
        threshold = int(sys.argv[3])
else:
    print "FORMAT: python AHE_taxatrim.py [option: -t (taxa), -l (length)] [folder with fasta] [trimlist or threshold]"
    print "EXAMPLE: python AHE_taxatrim.py -t ./fasta list.lst"
    print "EXAMPLE: python AHE_taxatrim.py -l ./fasta 100"
    sys.exit()

translist = []

if option == "-t":
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
            if option == "-t":
                if seq.id in translist:
                    locilist.add(infile)
                    break
            elif option == "-l":
                if len(seq.seq) > threshold:
                    locilist.add(infile)
                    break
            else:
                print "test"
        # if option == "-t":
        # 	if seq.id in translist:
        # 		locilist.add(infile)
        #         break
        # elif option == "-l":
        #     if len(seq.seq) > threshold:
        #         locilist.add(infile)
        #         break
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
        shutil.copy2(f+"/"+locusfname, "./subset")
print "done"
