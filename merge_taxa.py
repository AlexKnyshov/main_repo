from Bio import SeqIO
import glob
import os
import shutil
import sys
if len(sys.argv) == 3:
    opt = sys.argv[1]
    folder = sys.argv[2]
else:
    print "FORMAT: python merge_taxa.py [opt: -n (normal), -l (ARLemmon pipeline), -i (index first, for I/O issues)] [folder]"
    print "EXAMPLE: python merge_taxa.py -n ./fasta"
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

if opt == "-i":
    indices = {}
    locinames = set([])
    print "index files"
    for f in files:
        fnew = ".".join(f.split("/")[-1].split(".")[:-1])
        prog = "processing "+str(fnew)+"..."
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        indices[fnew] = SeqIO.index(f, "fasta")
        for seq in indices[fnew]:
            locinames.add(seq)
    print "indexed, Nloci:", len(locinames)
    print "write loci"
    for loc in locinames:
        prog = "processing "+loc
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        with open("./merged/"+loc+".fas", "w") as fhandle:
            for key, val in indices.items():
                if loc in val:
                    print >> fhandle, ">"+key+"\n"+str(indices[key][loc].seq).upper()

else:
    for f in files:
        fnew = f.split("/")[-1]
        locusfname = fnew.split("_")[-1].split(".")[0] #locus name
        align1 = SeqIO.parse(f, "fasta") # select seqs to copy
        prog = "processing "+str(fnew)+"..."
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        for seq in align1:
            if opt == "-l":
                fhandle = open("./merged/"+seq.id+".fasta", "a")
                print >> fhandle, ">"+locusfname+".1"
                print >> fhandle, str(seq.seq).replace("-", "").upper()
            else:
                fhandle = open("./merged/"+seq.id+".fas", "a")
                print >> fhandle, ">"+locusfname+"\n"+str(seq.seq ).replace("-", "").upper()
            fhandle.close()

print "done"