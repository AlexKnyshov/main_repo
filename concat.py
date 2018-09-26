from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
#from Bio.Alphabet.IUPAC import *
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
import csv
import sys
import glob
import os

def most_common(lst):
    return max(set(lst), key=lst.count)

def ss_parser(filename):
	output = {}
	with open(filename, "rU") as fcon:
		next(fcon)#skip xread
		next(fcon)#skip taxa / char counters
		for line in fcon:
			line=line.strip()
			if line == ";":
				break
			output[line.split(" ")[0]] = line.split(" ")[-1]
	return output

def check_alphabet(filename, fformat):#code for this function is modified from https://stackoverflow.com/questions/41588725/estimate-alphabet-in-biopython-from-fasta-file/41596563#41596563
	alphabets = [IUPAC.ambiguous_dna, IUPAC.protein]#[ambiguous_dna, extended_protein]#, unambiguous_dna, extended_dna, ambiguous_rna, unambiguous_rna]
	first_record = list(SeqIO.parse(filename, fformat))[0]
	detected = ""
	#check DNA first:
	leftover = set(str(first_record.seq).upper()) - set(alphabets[0].letters) - set(["-", "?"])
	if len(leftover) == 0:
		detected = "DNA"
	else:
		leftover = set(str(first_record.seq).upper()) - set(alphabets[1].letters) - set(["-", "?","X"])
		if len(leftover) == 0:
			detected = "Prot"
		else:
			print filename, "Error: unknown alphabet, problematic symbols:", leftover
			sys.exit()
	return detected

if len(sys.argv) == 3:
	inputfolder = sys.argv[1]
	partnum = sys.argv[2]

else:
	print "FORMAT: python concat.py [folder with fasta] [split to codon positions: -3 (yes), -1 (no)]"
	print "EXAMPLE: python concat.py ./fasta -1"
	print "output is written to COMBINED.phy, partitions are written to partitions.prt"
	sys.exit()

fasta = ["fasta", "fas", "fa"]
phylip = ["phylip", "phy"]
nexus = ["nexus", "nex"]

print "first pass: creating a list of taxa..."
d = {}
files = glob.glob(inputfolder+"/*")
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	if f.split(".")[-1] in fasta or f.split(".")[-1] in phylip or f.split(".")[-1] in nexus:
		if f.split(".")[-1] in fasta:
			fformat = "fasta"
		elif f.split(".")[-1] in nexus:
			fformat = "nexus"
		else:
			fformat = "phylip-relaxed"
		alph = check_alphabet(f, fformat)
		if alph == "DNA":
			alignment = AlignIO.read(f, fformat, alphabet = Gapped(IUPAC.ambiguous_dna))
		elif alph == "Prot":
			alignment = AlignIO.read(f, fformat, alphabet = Gapped(IUPAC.protein, '-'))
		length = alignment.get_alignment_length()
		for seq in alignment:
			if seq.id in d:
				d[seq.id].append(fn)
			else:
				d[seq.id] = []
				d[seq.id].append(fn)
	elif f.split(".")[-1]=="ss":
		matrix = ss_parser(f)
		for rec in matrix.keys():
			if rec in d:
				d[rec].append(fn)
			else:
				d[rec] = []
				d[rec].append(fn)
print len(d), "taxa found in input files"

print "second pass: concatenating..."
final_matrix = dict.fromkeys(d.keys(),"")
start = 0
end = 0
outputfile = open("partitions.prt", "w")
for f in files:
	fnew = f.split("/")
	fn = fnew[len(fnew)-1]
	length = 0
	missed = list(d)
	model = ""
	if f.split(".")[-1] in fasta or f.split(".")[-1] in phylip or f.split(".")[-1] in nexus:
		if f.split(".")[-1] in fasta:
			fformat = "fasta"
		elif f.split(".")[-1] in nexus:
			fformat = "nexus"
		else:
			fformat = "phylip-relaxed"
		alph = check_alphabet(f, fformat)
		if alph == "DNA":
			alignment = AlignIO.read(f, fformat, alphabet = Gapped(IUPAC.ambiguous_dna))
			model = "DNA"
		elif alph == "Prot":
			alignment = AlignIO.read(f, fformat, alphabet = Gapped(IUPAC.protein, '-'))
			model = "WAG"
		length = alignment.get_alignment_length()
		for seq in alignment:
			if seq.id in final_matrix:
			 	final_matrix[seq.id] += str(seq.seq)
			else:
				final_matrix[seq.id] = ""
				final_matrix[seq.id] += str(seq.seq)
			missed.remove(seq.id)
	elif f.split(".")[-1]=="ss":
		matrix = ss_parser(f)
		length = len(matrix.values()[0])
		model = "MULTI"
		for rec in matrix.keys():
			if rec in final_matrix:
				final_matrix[rec] += matrix[rec]
			else:
				final_matrix[rec] = ""
				final_matrix[rec] += matrix[rec]
			missed.remove(rec)
	else:#skip non fasta / ss files
		continue
	if len(missed) > 0:
		for m in missed:
			if m in final_matrix:
				final_matrix[m] += "?"*length
			else:
				final_matrix[m] = ""
				final_matrix[m] += "?"*length
	end = start + length - 1
	if partnum == "-3":
		if model == "DNA":
			print >> outputfile, model+", "+fn+"-1="+str(start+1)+"-"+str(end+1)+"\\3"
			print >> outputfile, model+", "+fn+"-2="+str(start+2)+"-"+str(end+1)+"\\3"
			print >> outputfile, model+", "+fn+"-3="+str(start+3)+"-"+str(end+1)+"\\3"
		else:
			print >> outputfile, model+", "+fn+"="+str(start+1)+"-"+str(end+1)
	elif partnum == "-1":
	 	print >> outputfile, model+", "+fn+"="+str(start+1)+"-"+str(end+1)
	prog = "working on partition "+str(fn)+": starts "+str(start+1)+", ends "+str(end+1)
  	sys.stdout.write(prog+"\r")
  	sys.stdout.flush()
	start = end + 1
outf = open("COMBINED.phy", "w")
print >> outf, str(len(final_matrix))+" "+str(start)
for rec in final_matrix.keys():
	print >> outf, str(rec)+" "+str(final_matrix[rec])
outf.close()
print "\ndone"
