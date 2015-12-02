import sys
import csv
filename = sys.argv[1]
values = []
loci = {}
with open(filename, "rb") as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		print row[0], row[1]
		loci[row[0]] = float(row[1])
		values.append(float(row[1]))
mn = min(values)
mx = max(values)
avg = sum(values)/len(values)
print mn, mx, avg
lower = float(raw_input("enter lower boundary:"))
upper = float(raw_input("enter upper boundary:"))
for x, y in loci.items():
	if lower < y < upper:
		print x