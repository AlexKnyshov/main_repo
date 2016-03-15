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
	ref_dict.add(refseq.id)

print "creating a new folder", outfolder
if not os.path.exists (outfolder):
    os.makedirs(outfolder) #creating folder if necessary

print "parsing mfs..."
for f in files:
	print f.split("/")[-1]
	input_handle = open(f, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	newal = []
	for seq in alignments:
		if seq.id in ref_dict:
			print seq.id, "found"
			newal.append(seq)
	outfile = open(outfolder+f.split("/")[-1], "w")
	SeqIO.write(newal, outfile, "fasta")
	outfile.close()

print "done"