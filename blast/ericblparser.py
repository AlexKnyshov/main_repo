from Bio import SeqIO
import glob
import os
import csv
import sys
blastfilearg = sys.argv[1]
trif = sys.argv[2]
evalue = float(sys.argv[3])

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
print "scanning the transcriptome..."
#output2 = []
warninglist = []
fhandle = open("result.fas", "w")
c1 = 0
for trifile in glob.glob(trif+"/*.fasta"):
    inputf = SeqIO.parse(trifile, "fasta")
    print "searching for contigs in:", trifile
    for seq in inputf:
        #print count
        if seq.id in output: #if contig is in ahe (was found as a blast hit)
            print "found", seq.id, "length:", len(seq.seq)
            SeqIO.write(seq, fhandle, "fasta")
            c1 += 1
fhandle.close()
print c1, "loci extracted", count, "original targets"
print "done"