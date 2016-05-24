from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment
import re
import sys
import glob
import os

if len(sys.argv) > 1:
	inputfolder = sys.argv[1]
else:
	print "FORMAT: python concat.py [folder with fasta] [split to codon positions: -3 (yes), -1 (no)] [partition_finder output: -pf2y, -pf2n] [phylip type: -i (interleaved), -s (sequential)]"
	print "EXAMPLE: python concat.py ./fasta -1 -i -pf2n"
	print "EXAMPLE: python concat.py ./fasta -1 -s -pf2y"
	print "output is written to COMBINED.phy, partitions are written to partitions.prt"
	sys.exit()

p1 = re.compile("---")
p2 = re.compile("-..")
p3 = re.compile("--.")
p4 = re.compile("--.")
p5 = re.compile(".-.")
p6 = re.compile(".--")
p7 = re.compile("..-")

goodcount = {}
badcount = []
print "input folder", inputfolder
files = glob.glob(inputfolder+"/*.fas")
if not os.path.exists ("./translated"):
	os.makedirs("./translated")
for f in files:
	goodlocus = False
	transseq = ""
	nuclseq = ""
	alignment = AlignIO.read(f, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
	for seq in alignment:
		t1 = str(seq.seq).replace("-", "N")
		t2 = t1.replace("?", "N")
		seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
		if seq.seq[0:3] != "NNN":
			nuclseq = seq.seq
			break
	if nuclseq != "":
		transseq = nuclseq.translate() #frame1
		print "locus", f
		if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 1
			#print "Nooo"
			transseq = nuclseq[1:].translate() #frame2
			if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 2
				transseq = nuclseq[2:].translate() #frame3
				if "*" in transseq[:-1] or "X" in transseq[:-1]: #check frame 3
					print "bad locus"
					badcount.append(f)
				else:
					goodlocus = True
					goodcount[f] = 3
			else:
				goodlocus = True
				goodcount[f] = 2
		else:
			goodlocus = True
			goodcount[f] = 1
	if goodlocus == True:
		goodlocus = False
		fnew = f.split("/")
		fn = fnew[len(fnew)-1]
		fn2 = "./translated/"+fn.split(".")[0]+".fas"
		outfile = open(fn2, "w")
		for seq in alignment:
			t1 = str(seq.seq).replace("-", "N")
			t2 = t1.replace("?", "N")
			seq.seq = Seq(t2, alphabet=Gapped(IUPAC.ambiguous_dna))
			if goodcount[f] == 1:
				transseq = seq.seq.translate()
			if goodcount[f] == 2:
				transseq = seq.seq[1:].translate()
			if goodcount[f] == 3:
				transseq = seq.seq[2:].translate()
			unambig_translate = str(transseq).replace("B", "X")
			unambig_translate = unambig_translate.replace("Z", "X")
			unambig_translate = unambig_translate.replace("J", "X")
			unambig_translate = unambig_translate.replace("U", "X")
			unambig_translate = unambig_translate.replace("O", "X")
			print >> outfile, ">"+seq.id, "\n", unambig_translate
		outfile.close()

print "good loci:", len(goodcount)
print "bad loci:"
for x in badcount:
	print x


print "done"
