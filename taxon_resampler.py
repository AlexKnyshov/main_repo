from Bio import SeqIO
from Bio import AlignIO
import sys
import glob
import os
import operator

print "init..."
folder = sys.argv[1] #folder with mfs
ext = sys.argv[2] #ext of mfs
files = glob.glob(folder+"/*"+ext)
infile = sys.argv[3] #old phylip file
outfolder = sys.argv[4] #new folder with trimmed taxon sets

ref_dict = set()
print "trying open phylip..."
ref_align = AlignIO.read(infile, "phylip-relaxed")
for refseq in ref_align:
	ref_dict.add(refseq.id.split("_")[-1])
print ref_dict

print "creating a new folder", outfolder
if not os.path.exists (outfolder):
    os.makedirs(outfolder) #creating folder if necessary

print "parsing mfs..."
l = []
for f in files:
	print f.split("/")[-1]
	input_handle = open(f, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	newal = []
	for seq in alignments:
		if seq.id.split("_")[-1] in ref_dict:
			#print seq.id, "found"
			newal.append(seq)
			l.append(seq.id.split("_")[-1])
	outfile = open(outfolder+f.split("/")[-1], "w")
	SeqIO.write(newal, outfile, "fasta")
	outfile.close()
	#print newal
l = set(l)
print l
print "checking reference against the extracted"
for item in ref_dict:
	if item not in l:
		print item, " was not extracted (-)"
	else:
		print item, " was extracted (+)"
print "done"