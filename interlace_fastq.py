from Bio import SeqIO
import sys


if len(sys.argv) == 3:
	readname1 = sys.argv[1]
	readname2 = sys.argv[2]

else:
	print "FORMAT: argument1 = read1, argument2 = read2"
	print "EXAMPLE: read1.fq read2.fq"
	print "output is written to reads.fq"
	sys.exit()

readhandle1 = open(readname1)
read1 = SeqIO.parse(readhandle1, "fastq")
# readhandle2 = open(readname2)
read2 = SeqIO.index(readname2, "fastq")
readhandle3 = open("reads.fq","w")
for seq in read1:
	SeqIO.write(seq, readhandle3, "fastq")
	SeqIO.write(read2[seq.id[:-2]+"/2"], readhandle3, "fastq")
readhandle1.close()
# readhandle2.close()
readhandle3.close()