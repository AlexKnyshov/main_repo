from Bio import Phylo
from Bio import SeqIO
import shutil
import os
import sys
import glob
import re

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

locusname = re.compile(".*(T.{1,3}_L.{1,4}fas).*")

if len(sys.argv) >= 3:
	files = glob.glob(sys.argv[2]+"/RAxML*")
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
	# print f
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
		treeavg = tree.total_branch_length()/len(tree.get_terminals())
		excludelist = []
		count = 0
		for clade in tree.find_clades():
			if len(tree.get_path(clade)) > 1:
				parent = get_parent(tree, clade)
				if tree.distance(clade, parent) / treeavg > 3:
					excludelist.append(clade)#.get_terminals()[0])
			else:
				if tree.distance(clade) / treeavg > 3:
					excludelist.append(clade)#.get_terminals()[0])

		if len(excludelist)>0:
			for item in excludelist:
				for term in item.get_terminals():
					trimlist.add(term.name)
			names = list(trimlist)
			result[f] = names
			#fname = locusname.match(f)
			fname = ".".join(f.split("/")[-1].split(".")[1:])
			aliout = open("./reduced/"+fname, "w")
			aliin = open(alifiles+"/"+fname, "rU")
			for seq in SeqIO.parse(aliin, "fasta"):
				if seq.id not in names:
				 	print >> aliout, ">"+seq.id, "\n", seq.seq
		 	status = "long branches detected"
			aliout.close()
			aliin.close()
		else:
			status = "no long branches detected"
	#progress bar
	progbarc +=1
	progbar = int(round(float(progbarc)/len(files)*100, 0))
	hashes = '#' * int(progbar * 0.2)
	spaces = ' ' * (20 - len(hashes))
	prog = "Progress: ["+str(hashes)+str(spaces)+"] "+str(progbar)+" %  "+str(status)+", total tree length: "+str(tree.total_branch_length())
 	sys.stdout.write(prog+"\r")
 	sys.stdout.flush()
#output the final table
with open("treesfilter.tab", "w") as outfile:
	for tc, tcv in result.items():
		print >> outfile, tc, "\t", tcv
print "Done"