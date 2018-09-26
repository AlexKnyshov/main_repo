from Bio import SeqIO
import glob
import os
import csv
import sys

if len(sys.argv) == 5:
    blastfilearg = sys.argv[1]
    assemblyf = sys.argv[2]
    evalue = float(sys.argv[3])
    opt = sys.argv[4]
else:
    print "FORMAT: python ericblparser.py [blast file] [folder or fasta file] [evalue] [option: -f (file), -d (folder)]"
    print "EXAMPLE: python ericblparser.py blast.out assembly.fasta 1e-10 -f"
    sys.exit()

output = set() #main dctionary
print "reading blastfile...", blastfilearg
blastfile = open(blastfilearg, "rU")
reader = csv.reader(blastfile, delimiter='\t')
for row in reader:
    if float(row[10]) <= evalue:
        output.add(row[1])
blastfile.close()
count = int(len(output))
print count, "targets found to be extracted"
print "extraction..."
fhandle = open("result.fas", "w")
c1 = 0
c2 = count
if opt == "-d":
    fls = glob.glob(assemblyf+"/*.fasta")
elif opt == "-f":
    fls = [assemblyf]
for trifile in fls:
    inputf = SeqIO.parse(trifile, "fasta")
    print "searching for contigs in:", trifile
    for seq in inputf:
        if seq.id in output:
            print "found", seq.id, "length:", len(seq.seq)
            SeqIO.write(seq, fhandle, "fasta")
            c1 += 1
            c2 -= 1
        if c2 == 0:
            break
    if c2 == 0:
        break
fhandle.close()
print c1, "contigs extracted", count, "original targets"
print "done"