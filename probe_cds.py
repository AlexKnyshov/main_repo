from Bio import SeqIO
import sys
import os
import shutil

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

def determine_direction(a,b,c,d):
	if a < b and c < d:
		return True
	elif a < b and c > d:
		return False
	elif a > b and c < d:
		return False
	elif a > b and c > d:
		return True

gffout = parse_gff(gfffile)

mkdirfunc("CDS")

with open(blastfile,"rU") as blasthandle:
	reader = csv.reader(blasthandle, delimiter='\t')
	for row in reader:
		contig = row[1]
		if contig in gffout:
			for cds in gffout[contig]:
				if getOverlap(cds[0:2],[int(row[8]), int(row[9])]) > 0:
					hitdirection = determine_direction(int(row[6]),int(row[7]),int(row[8]),int(row[9]))
					targets = min(int(row[8]), int(row[9]))
					targete = max(int(row[8]), int(row[9]))
					querys = min(int(row[6]), int(row[7]))-1
					querye = max(int(row[6]), int(row[7]))-1
					if cds[0] < targets:
						startgap = 0
					else:
						startgap = cds[0] - targets
					if cds[1] > targete:
						endgap = 0
					else:
						endgap = targete - cds[1]
					with open(row[0]) as infhandle:
						seq = SeqIO.read(infhandle, "fasta")
						if hitdirection:
							outseq = seq[(querys+startgap):(querye-endgap)].translate()
						else:
							outseq = seq[(querys+endgap):(querye-startgap)].reverse_complement().translate()
						basename = row[0].split("/")[-1]
						with open("./CDS/"+basename, "w") as outh:
							print >> outh, ">"+outseq.id
							print >> outh, outseq.seq