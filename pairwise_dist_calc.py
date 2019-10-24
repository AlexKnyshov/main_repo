from Bio import AlignIO
import sys
import glob

inputfolder = sys.argv[1]
files = glob.glob(inputfolder+"/*.fas")
list_pc_match = []
list_pc_mismatch = []
list_pc_bait_gap = []
list_pc_target_gap = []
for f in files:
	aln = AlignIO.read(open(f), "fasta")
	for seq in aln:
		seq.seq = seq.seq.upper()
	length = aln.get_alignment_length()
	bait_gap = 0
	target_gap = 0
	mismatch = 0
	match = 0
	for pos in range(length):
		if aln[0][pos] == "-" or aln[0][pos] == "n" or aln[0][pos] == "N":
			bait_gap += 1
		elif aln[1][pos] == "-" or aln[1][pos] == "n" or aln[1][pos] == "N":
			target_gap += 1
		else:
			if aln[0][pos] == aln[1][pos]:
				match += 1
			else:
				mismatch += 1
	pc_match = round(match / float(length) * 100, 3)
	list_pc_match.append(pc_match)
	pc_mismatch = round(mismatch / float(length) * 100, 3)
	list_pc_mismatch.append(pc_mismatch)
	pc_bait_gap = round(bait_gap / float(length) * 100, 3)
	list_pc_bait_gap.append(pc_bait_gap)
	pc_target_gap = round(target_gap / float(length) * 100, 3)
	list_pc_target_gap.append(pc_target_gap)
	print f+'\t'+str(pc_match)+'\t'+str(pc_mismatch)+'\t'+str(pc_bait_gap)+'\t'+str(pc_target_gap)

print "mean"+'\t'+str(sum(list_pc_match)/len(list_pc_match))+'\t'+str(sum(list_pc_mismatch)/len(list_pc_match))+'\t'+str(sum(list_pc_bait_gap)/len(list_pc_match))+'\t'+str(sum(list_pc_target_gap)/len(list_pc_match))
