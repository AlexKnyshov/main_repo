from Bio import Phylo
import os
import sys
import csv

if len(sys.argv) == 3:
	csvname = sys.argv[1]
	treename = sys.argv[2]

else:
	print "FORMAT: python relabelPhylo.py [names table] [tree file]"
	print "EXAMPLE: python relabelPhylo.py names.csv tree.tre"
	sys.exit()

renamer = {}
with open(csvname) as csvhandle:
	for row in csv.reader(csvhandle):
		renamer[row[0].strip()] = row[1].strip()

tree = Phylo.read(treename, "newick")
for tip in tree.get_terminals():
	if tip.name in renamer:
		tip.name = renamer[tip.name]
Phylo.write(tree, treename+".renamed", "newick")