from Bio import SeqIO
import glob
import os
import shutil
import sys
if len(sys.argv) == 3:
    folder1 = sys.argv[1]
    folder2 = sys.argv[2]
else:
    print "FORMAT: python merge_loci.py [folder1] [folder2]"
    print "EXAMPLE: python merge_loci.py ./fasta1 ./fasta2"
    sys.exit()

print "creating a list of taxa..."
files = {}
d = {}
files1 = glob.glob(folder1+"/*.fas")
files2 = glob.glob(folder2+"/*.fas")

if not os.path.exists ("./modified/"):
    os.makedirs("./modified") #creating folder if necessary
else:
    shutil.rmtree("./modified/") #removing old files
    os.makedirs("./modified")

#copy files
print "copying files:"
for x in files1:
    locusfname = x.split("/")[-1].split("_")[-1] #stripping off the experiment tag
    files[locusfname] = 0
    if not os.path.exists ("./modified/"+locusfname):
        prog = "copying "+str(locusfname)+"..."
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        shutil.copy2(x, "./modified/"+locusfname)
print ""

files1 = glob.glob("./modified/*.fas") #getting a new reference filelist

for f in files2: #checking folder2
    fnew = f.split("/")[-1] #filename of folder2
    locusfname = fnew.split("_")[-1] #locus name
    if "./modified/"+locusfname in files1: #if locus in reference folder
        align1 = SeqIO.parse("./modified/"+locusfname, "fasta") # select seqs to copy
        for ali1seq in align1:
            d[ali1seq.id] = 0
        align2 = SeqIO.parse(f, "fasta")
        fhandle = open("./modified/"+locusfname, "a") # copying the seqs
        for seq in align2:
            if seq.id in d:
                print "seq is present"
            else:
                print "writing seq..."
                SeqIO.write(seq, fhandle, "fasta")
        fhandle.close()
    else:
        print "copying entire locus..."
        shutil.copy2(folder2+"/"+fnew, "./modified/"+locusfname) #not in reference folder, copy

print "done"