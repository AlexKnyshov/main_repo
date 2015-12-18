from Bio import SeqIO
import csv
import sys
if len(sys.argv) == 4:
    blastfilearg = sys.argv[1]
    fastafilearg = sys.argv[2]
    outfilearg = sys.argv[3]
else:
    print "FORMAT: python blparser.py [blastfile] [fastafile] [ouputfasta]"
    print "EXAMPLE: python blparser.py results.out CespC_062414_SOAPtrans_K49.scafSeq test.fas"
    sys.exit()
blastfile = open(blastfilearg, "r")
output = {}
reader = csv.reader(blastfile, delimiter='\t')
currentkey = ""
currentmatch = ""
currente = 0.0
for row in reader:
    if currentkey != row[0]:
        if float(row[10]) <= 1e-40:
            if row[1] in output:
                print "warning: the key exists", row[1]
            output[row[1]] = row[10]
            print row[0], row[1], row[2], row[10], row[11]
    else:
        if currentmatch != row[1] and float(row[10]) <= 1e-40:
            print "warning: several matches detected", row[1], row[10], "delta is", currente-float(row[10])
    currentkey = row[0]
    currentmatch = row[1]
    currente = float(row[10])
blastfile.close()
#scanning the source
count = int(len(output))
output2 = []
inputf = SeqIO.parse(fastafilearg, "fasta")
print "searching for contigs in", fastafilearg
print "query length:", count

for seq in inputf:
    if seq.id in output:
        output2.append(seq)
        count -= 1
        print "found", seq.id
    elif count == 0:
        print "search terminated"
        break
SeqIO.write(output2, outfilearg, "fasta")