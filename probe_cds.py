from Bio import SeqIO
import sys
import os
import shutil
import csv

if len(sys.argv) >= 2:
	blastfile = sys.argv[1]
	gfffile = sys.argv[2]

else:
	print "FORMAT: argument1 = blast table, argument2 = gff3 file"
	print "EXAMPLE: ./blast.blast ./gff3.gff3"
	print "output is written to ./CDS folder"
	sys.exit()

def mkdirfunc(dir1):
	if not os.path.exists (dir1):
		os.makedirs(dir1) #creating folder if necessary
	else:
		shutil.rmtree(dir1) #removing old files
		os.makedirs(dir1)

def parse_gff(gffname):
	with open(gffname) as gffhandle:
		gffdict = {}
		for line in gffhandle:
			if line[0] != "#":
				line = line.strip().split()
				if line[2] == "CDS":
					contig = line[0]
					coord1 = int(line[3])
					coord2 = int(line[4])
					if line[6] == "+":
						strand = True
					else:
						strand = False
					cds = coord1, coord2, strand
					if contig in gffdict:
						gffdict[contig].append(cds)
					else:
						gffdict[contig] = [cds]
	return gffdict

def getOverlap(a, b):
	a0=min(a)
	a1=max(a)
	b0=min(b)
	b1=max(b)
	return max(0, min(a1, b1) - max(a0, b0))

def determine_direction(a,b):
	if a < b:
		return True
	else:
		return False

gffout = parse_gff(gfffile)

mkdirfunc("CDS")

badlist = {}
with open(blastfile,"rU") as blasthandle:
	reader = csv.reader(blasthandle, delimiter='\t')
	for row in reader:
		print row[0]
		contig = row[1]
		if contig in gffout:
			ovlp = False
			for cds in gffout[contig]:
				if getOverlap(cds[0:2],[int(row[8]), int(row[9])]) > 0:
					ovlp = True
					targets = min(int(row[8]), int(row[9]))
					targete = max(int(row[8]), int(row[9]))
					targetdirection = determine_direction(int(row[8]),int(row[9]))
					querys = min(int(row[6]), int(row[7]))-1
					querye = max(int(row[6]), int(row[7]))-1
					querydirection = determine_direction(int(row[6]),int(row[7]))
					if cds[0] < targets:
						if querydirection == targetdirection:
							startgap = 0
						else:
							endgap = 0
					else:
						if querydirection == targetdirection:
							startgap = cds[0] - targets
						else:
							endgap = cds[0] - targets
					if cds[1] > targete:
						if querydirection == targetdirection:
							endgap = 0
						else:
							startgap = 0
					else:
						if querydirection == targetdirection:
							endgap = targete - cds[1]
						else:
							startgap = targete - cds[1]
					cdsdirection = cds[2]
					# print cds, row[6:10]
					# print "ranges", querys, querye, targets, targete
					# print "startgap on target", startgap, "end", endgap
					# print "directions:", querydirection, targetdirection, cdsdirection
					with open(row[0]) as infhandle:
						seqs = list(SeqIO.parse(infhandle, "fasta"))
						seq = seqs[0].seq
						# print seq[(querys+startgap):(querye-endgap)].translate()
						if cdsdirection:
							if querydirection == targetdirection:
								outseq = seq[(querys+startgap):(querye-endgap)]
							else:
								outseq = seq[(querys+startgap):(querye-endgap)].reverse_complement()
						else:
							if querydirection == targetdirection:
								outseq = seq[(querys+startgap):(querye-endgap)].reverse_complement()
							else:
								outseq = seq[(querys+startgap):(querye-endgap)]
						basename = row[0].split("/")[-1]
						# if startgap > 0 or endgap > 0:
						# 	print basename, outseq, querys, querye
						with open("./CDS/"+basename, "w") as outh:
							print >> outh, ">"+seqs[0].id
							print >> outh, outseq
					break
			if not ovlp:
				badlist[row[0]] = "region not in CDS"
		else:
			badlist[row[0]] = "contig not in GFF"
	for key, val in badlist.items():
		print key, ":", val