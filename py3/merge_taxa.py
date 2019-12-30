from Bio import SeqIO
import glob
import os
import shutil
import sys
if len(sys.argv) == 3:
    opt = sys.argv[1]
    folder = sys.argv[2]
else:
    print ("FORMAT: python merge_taxa.py [opt: -n (normal), -l (ARLemmon pipeline)] [folder]")
    print ("EXAMPLE: python merge_taxa.py -n ./fasta")
    sys.exit()

if not os.path.exists ("./merged/"):
    os.makedirs("./merged") #creating folder if necessary
else:
    shutil.rmtree("./merged/") #removing old files
    os.makedirs("./merged/")

if opt  == "-l":
    files = glob.glob(folder+"/*.fasta")
else:
    files = glob.glob(folder+"/*.fas")

for f in files: #checking folder2
    fnew = f.split("/")[-1] #filename of folder2
    locusfname = fnew.split("_")[-1].split(".")[0] #locus name
    align1 = SeqIO.parse(f, "fasta") # select seqs to copy
    prog = "processing "+str(fnew)+"..."
    sys.stdout.write(prog+"\r")
    sys.stdout.flush()
    for seq in align1:
        if opt == "-l":
            fhandle = open("./merged/"+seq.id+".fasta", "a")
            print (">"+locusfname+".1", file=fhandle)
            print (str(seq.seq).replace("-", "").upper(), file=fhandle)
        else:
            fhandle = open("./merged/"+seq.id+".fas", "a")
            print (">"+locusfname, "\n", str(seq.seq ).replace("-", "").upper(), file=fhandle)
        fhandle.close()
print ("done")