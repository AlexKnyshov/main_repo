from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
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
			if "\t" in line:
				output[line.split("\t")[0]] = line.split("\t")[-1]
			else:
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

if len(sys.argv) >= 3:
	inputfolder = sys.argv[1]
	partnum = sys.argv[2]

else:
	print "FORMAT: python concat.py [folder with fasta] [split to codon positions: -3 (yes), -1 (no), -12 (combine first two), -12a (combine first two across all), -tnt, -nex, -nex2 (recode DNA as discrete)]"
	print "EXAMPLE: python concat.py ./fasta -1"
	print "output is written to COMBINED.phy, partitions are written to partitions.prt"
	sys.exit()

fasta = ["fasta", "fas", "fa"]
phylip = ["phylip", "phy"]
nexus = ["nexus", "nex"]

files = sorted(glob.glob(inputfolder+"/*"))

print "concatenating..."
final_matrix = {}
d = set([])
start = 0
end = 0
if partnum != "-tnt" and partnum != "-nex" and partnum != "-nex2":
	outputfile = open("partitions.prt", "w")
if partnum == "-12a":
	range12 = []
	range3 = []
models = []
starts = []
ends = []
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
			model = "LG"
		length = alignment.get_alignment_length()
		for seq in alignment:
			if seq.id in final_matrix:
			 	final_matrix[seq.id] += str(seq.seq)
			else:
				if len(ends) > 0:
					final_matrix[seq.id] = "?"*ends[-1]
				else:
					final_matrix[seq.id] = ""
				final_matrix[seq.id] += str(seq.seq)
				d.add(seq.id)
			if seq.id in missed:
				missed.remove(seq.id)
	elif f.split(".")[-1]=="ss":
		matrix = ss_parser(f)
		length = len(matrix.values()[0])
		model = "MULTI"
		for rec in matrix.keys():
			if rec in final_matrix:
				final_matrix[rec] += matrix[rec]
			else:
				if len(ends) > 0:
					final_matrix[rec] = "?"*ends[-1]
				else:
					final_matrix[rec] = ""
				final_matrix[rec] += matrix[rec]
				d.add(rec)
			if rec in missed:
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
	models.append(model)
	starts.append(start+1)
	ends.append(end+1)
	if partnum == "-3":
		if model == "DNA":
			print >> outputfile, model+", "+fn+"-1="+str(start+1)+"-"+str(end+1)+"\\3"
			print >> outputfile, model+", "+fn+"-2="+str(start+2)+"-"+str(end+1)+"\\3"
			print >> outputfile, model+", "+fn+"-3="+str(start+3)+"-"+str(end+1)+"\\3"
		else:
			print >> outputfile, model+", "+fn+"="+str(start+1)+"-"+str(end+1)
	elif partnum == "-12":
		if model == "DNA":
			print >> outputfile, model+", "+fn+"-12="+str(start+1)+"-"+str(end+1)+"\\3, "+str(start+2)+"-"+str(end+1)+"\\3"
			print >> outputfile, model+", "+fn+"-3="+str(start+3)+"-"+str(end+1)+"\\3"
		else:
			print >> outputfile, model+", "+fn+"="+str(start+1)+"-"+str(end+1)
	elif partnum == "-12a":
		if model == "DNA":
			range12.append(str(start+1)+"-"+str(end+1)+"\\3")
			range12.append(str(start+2)+"-"+str(end+1)+"\\3")
			range3.append(str(start+3)+"-"+str(end+1)+"\\3")
		else:
			print "model", model, "is not supported, opt -12a is only for DNA"
			sys.exit()
	elif partnum == "-1":
	 	print >> outputfile, model+", "+fn+"="+str(start+1)+"-"+str(end+1)
	prog = "working on partition "+str(len(starts))+", "+str(fn)+": starts "+str(start+1)+", ends "+str(end+1)
  	sys.stdout.write(prog+"\r")
  	sys.stdout.flush()
	start = end + 1
if partnum == "-tnt":
	outf = open("COMBINED.tnt", "w")
	print >> outf, "xread"
	print >> outf, str(start)+" "+str(len(final_matrix))
	tntm = []
	tnts = []
	tnte = []
	for part in range(0, len(models)):
		if part == 0:
			tnts.append(starts[0]-1)
			tnte.append(ends[0])
			tntm.append(models[0])
		else:
			if models[part] == models[part-1]:
				tnte[-1] = ends[part]
			else:
				tnts.append(starts[part]-1)
				tnte.append(ends[part])
				tntm.append(models[part])
	for tntp in range(0, len(tntm)):
		if tntm[tntp] == "DNA":
			print >> outf, "&[dna]"
		elif tntm[tntp] == "MULTI":
			print >> outf, "&[num]"
		for rec in sorted(final_matrix.keys()):
			print >> outf, str(rec)+" "+str(final_matrix[rec][tnts[tntp]:tnte[tntp]])
	print >> outf, ";"
	print >> outf, "proc/;"
	outf.close()
elif partnum == "-nex":
	outf = open("COMBINED.nex", "w")
	mline = ""
	for part in range(0, len(models)):
		if models[part] == "MULTI":
			models[part] = "Standard"
		if part == 0:
			mline += models[0]
			mline += ":"
			mline += str(starts[0])
		else:
			if models[part] != models[part-1]:
				mline += "-"
				mline += str(ends[part-1])
				mline += ","
				mline += models[part]
				mline += ":"
				mline += str(starts[part])
	mline += "-"
	mline += str(ends[part])
	print >> outf, "#nexus"
	print >> outf, "begin data;"
	print >> outf, "dimensions ntax="+str(len(final_matrix))+" nchar="+str(start)+";"
	print >> outf, "format datatype=mixed("+mline+") interleave=yes  GAP = - MISSING = ?;"
	print >> outf, "matrix"
	for rec in sorted(final_matrix.keys()):
		print >> outf, str(rec)+"\t"+str(final_matrix[rec])
	print >> outf, ";"
	print >> outf, "end;"
	for part in range(0, len(models)):
		print >> outf, "charset", models[part]+str(part)+"="+str(starts[part])+"-"+str(ends[part])+";"
	outf.close()
elif partnum == "-nex2":
	outf = open("COMBINED.nex", "w")
	print >> outf, "#nexus"
	print >> outf, "begin data;"
	print >> outf, "dimensions ntax="+str(len(final_matrix))+" nchar="+str(start)+";"
	print >> outf, "format datatype=STANDARD interleave=yes  GAP = - MISSING = ? SYMBOLS = \"  0 1 2 3 4 5\";"
	print >> outf, "matrix"
	for rec in sorted(final_matrix.keys()):
		print >> outf, str(rec)+"\t"+str(final_matrix[rec]).replace("a", "0").replace("t", "1").replace("g", "2").replace("c", "3").replace("w", "?").replace("n", "?").replace("r", "?").replace("y", "?").replace("s", "?").replace("m", "?")
	print >> outf, ";"
	print >> outf, "end;"
	outf.close()
else:
	outf = open("COMBINED.phy", "w")
	print >> outf, str(len(final_matrix))+" "+str(start)
	for rec in sorted(final_matrix.keys()):
		print >> outf, str(rec)+" "+str(final_matrix[rec])
	outf.close()
	if partnum == "-12a":
		print >> outputfile, model+", concat-12="+",".join(range12)
		print >> outputfile, model+", concat-3="+",".join(range3)
	outputfile.close()
print "\ndone"
