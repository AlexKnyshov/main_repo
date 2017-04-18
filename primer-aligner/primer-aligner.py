from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
import os
import csv
import sys
import re

def revcomfunc(seq):
	revcom = ""
	reverse = seq.upper()[::-1]
	for y in reverse:
		if y == "A":
			revcom+="T"
		elif y == "T":
			revcom+="A"
		elif y == "G":
			revcom+="C"
		elif y == "C":
			revcom+="G"
	return revcom

if len(sys.argv) == 5:
		subjectfile = sys.argv[1]
		primerfile = sys.argv[2]
		blastfile = sys.argv[3]
		collapseopt = sys.argv[4]
		if collapseopt == "-c":
			collapse = True
		elif collapseopt == "-n":
			collapse = False
		else:
			print "unrecognized option"
			sys.exit()
else:
	print "FORMAT: python primer-aligner.py [subject sequence in fasta] [fasta file with primers] [blast output file] [collapse: -c (yes), -n (no)]"
	print "EXAMPLE: python primer-aligner.py rRNA_operon.fas primers.fas blast.blast -n"
	sys.exit()
blast = open(blastfile, "rU")
reader = csv.reader(blast, delimiter='\t')
primerhandle = open(primerfile, "r")
primers = SeqIO.parse(primerhandle, "fasta")
subjecthandle = open(subjectfile, "r")
subject = SeqIO.read(subjecthandle, "fasta")
ugapsubject = subject.seq.ungap("-")
print "length of the subject sequence", len(subject)
primerdict = {}
for p in primers:
	primerdict[p.id] = str(p.seq)
outhandle = open("output-primer-align.fas", "w")
if collapse:
	megasec = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna))
else:
	print >> outhandle, ">"+subject.id, "\n", subject.seq
for row in reader:
	tempseq = ""
	print row[0], row[8], row[9]
	if int(row[8])>int(row[9]):
		q = re.search(str(ugapsubject)[int(row[9])-1:int(row[8])], str(subject.seq))
		new9 = q.start()-int(row[6])+1
		new8 = q.end()
		revcom = revcomfunc(primerdict[row[0]])
		primerseqct = 0
		print "reverse"
		for x in range(len(subject)):
			if x < new9:
				tempseq += "-"
			elif new9<=x<new8 and primerseqct < len(primerdict[row[0]]):
				tempseq += revcom[primerseqct]
				primerseqct += 1
			else:
				tempseq += "-"
	else:
		q = re.search(str(ugapsubject)[int(row[8])-1:int(row[9])], str(subject.seq))
		new8 = q.start()-int(row[6])+1
		new9 = q.end()
		primerseqct = 0
		print "forward"
		for x in range(len(subject)):
			if x < new8:
				tempseq += "-"
			elif new8<=x<new9 and primerseqct < len(primerdict[row[0]]):
				tempseq += primerdict[row[0]][primerseqct]
				primerseqct += 1
			else:
				tempseq += "-"
	if collapse:
		megasec.append(SeqRecord(Seq(tempseq, Gapped(IUPAC.ambiguous_dna)), id=row[0]))
	else:
		print >> outhandle, ">"+row[0], "\n", tempseq
if collapse:
	conseq = ""
	for x in range(megasec.get_alignment_length()):
		pos = megasec[:,x]
		c = len(pos)
		for y in pos:
			if y != "-":
				conseq += y
				break
			if c == 1:
				conseq += "-"
				break
			c -= 1
	print >> outhandle, ">"+subject.id, "\n", subject.seq
	print >> outhandle, ">primers", "\n", conseq
blast.close()
primerhandle.close()
subjecthandle.close()
outhandle.close()
print "done"