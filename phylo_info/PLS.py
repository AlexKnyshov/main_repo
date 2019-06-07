import sys

tab = sys.argv[1]
prts = sys.argv[2]

tabout = {}
# with open(tab, "rb") as tabfile:
# 	next(tabfile)
# 	for row in tabfile:
# 		row1 = row.strip().split("\t")
# 		tabout[int(row1[0][2:])] = [float(i) for i in row1[1].split(" ")]
# #print tabout["tr1"][1:4]

with open(tab, "rb") as tabfile:
	next(tabfile)
	for row in tabfile:
		row1 = row.strip().split()
		#print row1
		#break
		tabout[int(row1[0][4:])] = [float(i) for i in row1[1:]]
# print tabout[1][1:4]


prtout = {}
# with open(prts, "rb") as prtfile:
# 	for row in prtfile:
# 		row1 = row.strip().split(" ")
# 		prtout[row1[1]] = [int(i) for i in row1[3].split("-")]
# print prtout["L74.fas"][0]

with open(prts, "rb") as prtfile:
	for row in prtfile:
		row1 = row.strip().split(" ")[1].split("=")
		prtout[row1[0]] = [int(i) for i in row1[1].split("-")]
#print prtout["L74.fas"]


finaout = open("prtlls.csv","w")
#print >> finaout, "partition,"+",".join(sorted(tabout.keys()))
print >> finaout, "partition,"+",".join([str(i) for i in sorted(tabout.keys())])
for partition, prtrange in sorted(prtout.items()):
	treeprtll = []
	#print partition, prtrange[0]-1, prtrange[1]-1
	for tree, psll in sorted(tabout.items()):
		prtll = sum(psll[prtrange[0]-1:prtrange[1]-1])
		treeprtll.append(str(prtll))
	print >> finaout, partition+","+",".join(treeprtll)
finaout.close()
print "done"