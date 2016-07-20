from Bio import Phylo
from Bio import SeqIO
import shutil
import os
import sys
import glob
import re

locusname = re.compile(".*(T.{1,3}_L.{1,4}fas).*")

if len(sys.argv) >= 3:
	files = glob.glob(sys.argv[2]+"/*")
	opt = sys.argv[1]
	if opt == "-lb":
		alifiles = sys.argv[3]
		print "creating an output folder..."
		if not os.path.exists ("./reduced/"):
		    os.makedirs("./reduced") #creating folder if necessary
		else:
		    shutil.rmtree("./reduced/") #removing old files
		    os.makedirs("./reduced")

else:
	print "FORMAT: python treeprocsep.py [option: -75, -avg, -lb] [folder with trees] ([folder with alignments])"
	print "EXAMPLE: python treeprocsep.py -75 ./trees"
	sys.exit()
if len(files) == 0:
	print "no trees in the directory"
#starting to process files
progbarc = 0
result = {}
#test
for f in files:
	print f
	tree = Phylo.read(f, "newick")
	#tree = Phylo.read("./../testtree.tre", "newick")
	if opt == "-75" or opt == "-avg":

		if tree.is_bifurcating():
			print "bifurcating"
		else:
			print "not bifurcating"
		counter = 0
		ct = 0
		test = tree.as_phyloxml()
		list4 = test.find_clades()
		totconf = 0
		lowconf = 0
		conflist = []
		for l1 in list4:
			counter +=1
			if l1.is_bifurcating():
				ct +=1
			l2 = str(l1._get_confidence())
			if l2.find("Confidence") == 0:
				totconf +=1
				l3=l2.split(",")
				l4 = l3[1][:-1].split("=")
				if opt == "-75":
					if int(l4[1]) <75:
						lowconf +=1
				elif opt == "-avg":
					conflist.append(int(l4[1]))
		print "below 70 per cent", float(lowconf)/(totconf)
		print "clades #", counter, ct
		fname = f.split("/")[2]
		tr = fname.split(".")[1]
		if opt == "-75":
			result[tr] = float(lowconf)/(totconf)
		elif opt == "-avg":
			result[tr] = float(sum(conflist)/len(conflist))
	elif opt == "-lb":
		trimlist = set()
		print "total tree length", tree.total_branch_length()
		#result[f] = tree.total_branch_length()
		#terms = tree.get_terminals()
		#print terminal_neighbor_dists(tree)
		#for term in terms:
		#	print trimlist
		#	tpar = tree.collapse(term)
			#print tree.depths()
			#print tpar, tree.distance(tpar)
		# for x in tree.find_elements(terminal=True):
		# 	if tree.distance(x) / tree.total_branch_length() > 0.7:
		# 		result[f] = x.name
		# 		trimlist.add(x.name)
		excludelist = []
		for clade in tree.find_clades():
			if tree.distance(clade) / tree.total_branch_length() > 0.7:
				#print clade#.get_terminals()[0]
				excludelist.append(clade)#.get_terminals()[0])
				#result[f] = clade.get_terminals()
				#trimlist.add(clade.get_terminals())
		if len(excludelist)>0:
			print excludelist[0].get_terminals()
			names = []
			for name in excludelist[0].get_terminals():
				names.append(name.name)
			result[f] = names
			#print term, tree.distance(tree.find_clades(term, terminal=True))
			#print trimlist
		#for z in trimlist:
		#	print z
		# if len(trimlist)>0:
		# 	print trimlist
			fname = locusname.match(f)
			aliout = open("./reduced/"+fname.group(1), "w")
			#print trimlist
			aliin = open(alifiles+"/"+fname.group(1), "rU")
			#alignments = SeqIO.parse(aliin, "fasta")
			#print trimlist
			#if f in result:
			#	print result[f], "1"
			#trimlist.
			for seq in SeqIO.parse(aliin, "fasta"):
			#for seq in alignments:
				#print >> aliout, trimlist
				# if f in result:
				# 	print result[f], "2"
				if seq.id not in names:
					#SeqIO.write(seq, aliout, "fasta")
				 	# print seq.id, "NO"#, trimlist
				 	print >> aliout, ">"+seq.id, "\n", seq.seq
				else:
				 	print seq.id, "YES"
				 	#print >> aliout, seq.id, "YES"
			aliout.close()
			aliin.close()
		else:
			print "GOOD LOCUS"
	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	hashes = '#' * int(progbar * 0.2)
	spaces = ' ' * (20 - len(hashes))
	print "-------------------------------------"
	print "\rProgress: [{0}] {1}%".format(hashes + spaces, progbar)
	print "-------------------------------------"
#output the final table
with open("treesfilter.tab", "w") as outfile:
	for tc, tcv in result.items():
		print >> outfile, tc, "\t", tcv
print "Done"