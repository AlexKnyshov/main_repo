from Bio import SeqIO
#print "test"
infile = open("./fasta/T58_L1.phylip.fas", "r")
seqs = {}
for seq in SeqIO.parse(infile, "fasta"):
    seqs[seq.id] = str(seq.seq)
#pairwise
names =[]
for key in seqs.keys():
	names.append(key)
distlist = []
for num in range(len(names)-1):
	#print num
	pair1 = []
	pair2 = []
	pair1 += seqs[names[num]]
	pair2 += seqs[names[num+1]]
	#print pair1
	print pair2
	counter = 0
	countert = 0
	counterindel = 0
	startpos = 0
	endpos = len(pair1)
	#forward trim
	for index in range(len(pair1)):
		if pair1[index] != "-" or pair2[index] != "-":
			startpos = index
			break
	#reverse trim
	for index in range(len(pair1)-1, -1, -1):
		if pair1[index] != "-" or pair2[index] != "-":
			endpos = index
			break
	#calculating the distance
	for index in range(startpos, endpos+1):
		if pair1[index] == "A" and pair2[index] == "G":
			countert += 1
		elif pair1[index] == "G" and pair2[index] == "A":
			countert += 1
		elif pair1[index] == "A" and pair2[index] == "C":
			counter +=1
		elif pair1[index] == "A" and pair2[index] == "T":
			counter +=1
		elif pair1[index] == "T" and pair2[index] == "A":
			counter +=1
		elif pair1[index] == "T" and pair2[index] == "G":
			counter +=1
		elif pair1[index] == "C" and pair2[index] == "T":
			countert +=1
		elif pair1[index] == "T" and pair2[index] == "C":
			countert +=1
		elif pair1[index] == "C" and pair2[index] == "A":
			counter +=1
		elif pair1[index] == "C" and pair2[index] == "G":
			counter +=1
		elif pair1[index] == "G" and pair2[index] == "C":
			counter +=1
		elif pair1[index] == "G" and pair2[index] == "T":
			counter +=1
		elif pair1[index] == "-" or pair2[index] == "-":
			counterindel +=1
	#print float(counter)/(endpos-startpos), counter, countert, counterindel
	#check subtraction
	distlist.append(float(counter)/(endpos-startpos))
print sum(distlist)/len(distlist)
