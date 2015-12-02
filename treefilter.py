import sys
import csv

def hist(lst, bins=10):
	step = (mx-mn)/bins
	b = []
	first = True
	cmn = 0
	cmx = 0
	for x in range(bins):
		current = 0
		if first:
			for val in values:
				cmn = mn
				cmx = mn+step
				if cmn <= val < cmx:
					current +=1
			first = False
			cmn = cmn+step
			cmx = cmx+step
		else:
			for val in values:
				if cmn <= val < cmx:
					current +=1
			cmn = cmn+step
			cmx = cmx+step
		b.append(current)
	countbin = 1
	for bin in b:
		print round(mn+(step*countbin),2), "\t", "|", "#"*int(bin/float(sum(b))*100)
		countbin+=1

filename = sys.argv[1]
values = []
loci = {}
with open(filename, "rb") as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		#print row[0], row[1]
		loci[row[0]] = float(row[1])
		values.append(float(row[1]))
mn = min(values)
mx = max(values)
avg = sum(values)/len(values)
print "minimum", mn
print "maximum", mx
print "average", avg

hist(values,10)

lower = float(raw_input("enter lower boundary:"))
upper = float(raw_input("enter upper boundary:"))
outfile = open("filteredloci.tab", "w")
for x, y in loci.items():
	if lower < y < upper:
		print >> outfile, x+".phylip"
outfile.close()