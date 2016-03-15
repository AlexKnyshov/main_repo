from Bio import SeqIO
import glob
import os
import shutil
import csv
import sys
if len(sys.argv) == 4:
    blastfilearg = sys.argv[1]
    trif = sys.argv[2]
    ahefoldarg = sys.argv[3]
else:
    print "FORMAT: python blparser.py [blastfile] [asemblyfile] [ahefolder]"
    print "EXAMPLE: python blparser.py blast.tab trinity.fas ./fasta"
    sys.exit()

output = {} #main dctionary

#reading the blastfile
print "reading blastfile...", blastfilearg
blastfile = open(blastfilearg, "rU")
reader = csv.reader(blastfile, delimiter='\t')
currentkey = ""
for row in reader:
    for row in reader:
        if currentkey != row[0]: ##new query
            if float(row[10]) <= 1e-40:
                if row[0].split("//")[-1] in output: ##query is present
                    print "warning: the key exists", row[0].split("//")[-1]
                else:
                    output[row[0].split("//")[-1]] = row[1] ##standart output
                    print row[0].split("//")[-1], row[1], row[2], row[10], row[11]
        else: ##same query
            if currentmatch != row[1] and float(row[10]) <= 1e-40:
                print "warning: several matches detected", row[0].split("//")[-1].split(".tx_tm")[0], row[1], "delta is", currente-float(row[10])
        currentkey = row[0]
        currentmatch = row[1]
        currente = float(row[10])
blastfile.close()
count = int(len(output))
print count, "targets found to be extracted"

print output

#scanning the transcriptomes
if not os.path.exists ("./modified/"):
    os.makedirs("./modified") #creating folder if necessary
else:
    shutil.rmtree("./modified/") #removing old files
    os.makedirs("./modified")

print "scanning the transcriptome..."
#output2 = []
inputf = SeqIO.parse(trif, "fasta")
print "searching for contigs in:", trif
c1 = 0
for seq in inputf:
    if seq.id in output.values(): #if contig is in ahe
        print "start"
        for x,y in output.items(): #checking all ahe
            if seq.id == y:
                locusfname = x
                #print locusfname
                # if not os.path.exists ("./modified/"+locusfname):
                #     shutil.copy2(ahefoldarg+locusfname, "./modified")
                print "found", locusfname, y
                fhandle = open("./modified/"+locusfname, "a")
                seq.id = "ceratotrans"
                seq.name =""
                seq.description =""
                #SeqIO.write(seq, fhandle, "fasta")
                c1 += 1
            else:
                print x, y
        count -= 1
    elif count == 0:
        print "search terminated"
        break
print c1, "loci extracted"
count = int(len(output))

print "done"