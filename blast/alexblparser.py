from Bio import SeqIO
import glob
import os
import shutil
import csv
import sys
if len(sys.argv) == 4:
    blastfilearg = sys.argv[1]
    trif = sys.argv[2]
    #triextarg = sys.argv[3]
    # clcfoldarg = sys.argv[4]
    # clcextarg = sys.argv[5]
    ahefoldarg = sys.argv[3]
else:
    print "FORMAT: python blparser.py [blastfile] [trinityfolder] [triextarg] [clcfolder] [clcextarg] [ahefolder]"
    print "EXAMPLE: python blparser.py blast ./trinity fasta ./clc fa ./fasta"
    sys.exit()

output = {} #main dctionary
lociset = [] #list to count unique ahe loci processed

#reading the blastfile
print "reading blastfile", blastfilearg
blastfile = open(blastfilearg, "rU")
reader = csv.reader(blastfile, delimiter='\t')
currentkey = ""
for row in reader:
    for row in reader:
        if currentkey != row[0]:
            if float(row[10]) <= 1e-40:
                if row[1] in output:
                    print "warning: the key exists", row[1]
                output[row[0].split("//")[-1]] = row[1]
                print row[0].split("//")[-1], row[1], row[2], row[10], row[11]
        else:
            if currentmatch != row[1] and float(row[10]) <= 1e-40:
                print "warning: several matches detected", row[0].split("//")[-1].split(".tx_tm")[0], row[1], "delta is", currente-float(row[10])
        currentkey = row[0]
        currentmatch = row[1]
        currente = float(row[10])
        lociset.append(row[0])
    # if currentkey != row[0]+row[1][:5]:
    #     #if row[1] in output.values():
    #     for x, y in output.items():
    #         if y == row[1]:
    #             print "warning: the target exists:", y, "; the problem AHE locus is", x
    #     output[row[0]+row[1][:5]] = row[1]
    #     #print row[0], row[1], row[2], row[10], row[11]
    #     print row[0]+row[1][:5], row[1]#, row[1], currentkey
    # currentkey = row[0]+row[1][:5]
    # lociset.append(row[0]+row[1][:5])
blastfile.close()
count = int(len(output))
print count, "targets found to be extracted"


#scanning the transcriptomes
if not os.path.exists ("./modified/"):
    os.makedirs("./modified") #creating folder if necessary
else:
    shutil.rmtree("./modified/") #removing old files
    os.makedirs("./modified")

print "scanning the transcriptome..."
output2 = []
inputf = SeqIO.parse(trif, "fasta")
print "searching for contigs in", trif
c1 = 0
for seq in inputf:
    if seq.id in output.values():
        for x,y in output.items():
            if y == seq.id:
                locusfname = x
                print locusfname
                if not os.path.exists ("./modified/"+locusfname):
                    shutil.copy2(ahefoldarg+locusfname, "./modified")
                print "found", locusfname, y
                fhandle = open("./modified/"+locusfname, "a")
                seq.id = "ceratotrans"#seq.id[:5]
                seq.name =""
                seq.description =""
                SeqIO.write(seq, fhandle, "fasta")
        output2.append(seq)
        c1 += 1
        count -= 1
    elif count == 0:
        print "search terminated"
        break
print c1, "loci extracted"
count = int(len(output))

# print "scanning clc assemblies..."
# files = glob.glob(clcfoldarg+"/*."+clcextarg)
# for clcf in files:
#     output2 = []
#     inputf = SeqIO.parse(clcf, "fasta")
#     print "searching for contigs in", clcf
#     c1 = 0
#     for seq in inputf:
#         if seq.id in output.values():
#             for x,y in output.items():
#                 if y == seq.id:
#                     locusfname = x.split("fas")[0]+"fas"
#                     if not os.path.exists ("./modified/"+locusfname):
#                         shutil.copy2(ahefoldarg+locusfname, "./modified")
#                     print "found", locusfname, y
#                     fhandle = open("./modified/"+locusfname, "a")
#                     seq.id = seq.id[:5]
#                     seq.name =""
#                     seq.description =""
#                     SeqIO.write(seq, fhandle, "fasta")
#             output2.append(seq)
#             c1 += 1
#             count -= 1
#         elif count == 0:
#             print "search terminated"
#             break
#     print c1, "loci extracted"

print len(lociset), "total records processed"
print len(set(lociset)), "AHE loci found in transcriptomes"
print "done"