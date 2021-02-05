import sys
from Bio import AlignIO
from collections import Counter

infile = sys.argv[1]

endlen = 5
window = 30
winstep = 15
wingapthresh = 0.2 #amount of gaps to start checking seq for windows
gapthresh = 0.25 #amount of gaps within window to start detect ends 
threshold = 0.25 #freq below which base is outlier


trimdict = {}

with open(infile) as inhandle:
	alignment = AlignIO.read(inhandle, "fasta")
	allen = len(alignment[0,:])
	#bite global ends
	fineseqs = set()
	for p in range(endlen):
		posdict = Counter(alignment[:,p])
		totalnongap = float(len(alignment[:,p].replace("-","")))
		for s1 in range(len(alignment[:,p])):
			if totalnongap > 0:
				s = alignment[s1,p]
				if s != "-" and s!= "X":
					if alignment[s1].id not in fineseqs and posdict[s] / totalnongap < threshold:
						if alignment[s1].id in trimdict:
							trimdict[alignment[s1].id].add(p)
						else:
							trimdict[alignment[s1].id] = set([p])
					else:
						fineseqs.add(alignment[s1].id)
			else:
				if alignment[s1].id in trimdict:
					trimdict[alignment[s1].id].add(p)
				else:
					trimdict[alignment[s1].id] = set([p])
	fineseqs = set()
	for p in range(allen-1, allen-endlen-1, -1):
		posdict = Counter(alignment[:,p])
		totalnongap = float(len(alignment[:,p].replace("-","").replace("X","")))
		for s1 in range(len(alignment[:,p])):
			if totalnongap > 0:
				s = alignment[s1,p]
				if s != "-" and s!= "X":
					if alignment[s1].id not in fineseqs and posdict[s] / totalnongap < threshold:
						if alignment[s1].id in trimdict:
							trimdict[alignment[s1].id].add(p)
						else:
							trimdict[alignment[s1].id] = set([p])
					else:
						fineseqs.add(alignment[s1].id)
			else:
				if alignment[s1].id in trimdict:
					trimdict[alignment[s1].id].add(p)
				else:
					trimdict[alignment[s1].id] = set([p])

	#bite internal ends
	for s in alignment:
		seqdict = Counter(s)
		if seqdict["-"] / float(allen) + seqdict["X"] / float(allen) > wingapthresh:
			junctions = {}
			winbreak = False
			for w in range(0, allen, winstep):
				if w+window < allen:
					win = s[w:w+window]
				else:
					win = s[w:allen]
					winbreak = True
				windict = Counter(win)
				totalgap = 0
				if "-" in windict:
					totalgap += windict["-"]
				if "X" in windict:
					totalgap += windict["X"]
				if totalgap < len(win) and totalgap / float(len(win)) > gapthresh:
					for p in range(len(win)-1):
						p1 = win[p]
						p2 = win[p+1]
						if (p1 == "-" or p1 == "X") and (p2 != "-" and p2 != "X"):
							junctions[str(p+1+w)] = "F"
						elif (p1 != "-" and p1 != "X") and (p2 == "-" or p2 == "X"):
							junctions[str(p+w)] = "R"
				if winbreak:
					break
			for j in junctions:
				if junctions[j] == "F":
					if int(j)+endlen < allen:
						tailcoord = int(j)+endlen
					else:
						tailcoord = allen
					for p in range(int(j), tailcoord):
						posdict = Counter(alignment[:,p])
						totalnongap = float(len(alignment[:,p].replace("-","").replace("X","")))
						if totalnongap > 0:
							if posdict[s[p]] / totalnongap < threshold:
								if s.id in trimdict:
									trimdict[s.id].add((int(j),int(j)+endlen))
								else:
									trimdict[s.id] = set([(int(j),int(j)+endlen)])
								break
				if junctions[j] == "R":
					for p in range(int(j), int(j)-endlen, -1):
						posdict = Counter(alignment[:,p])
						totalnongap = float(len(alignment[:,p].replace("-","").replace("X","")))
						if totalnongap > 0:
							if posdict[s[p]] / totalnongap < threshold:
								if s.id in trimdict:
									trimdict[s.id].add((int(j)-endlen+1,int(j)+1))
								else:
									trimdict[s.id] = set([(int(j)-endlen+1,int(j)+1)])
								break

	#output trimmed
	trimposcounter = 0
	with open(infile+".masked", "w") as outhandle:
		for s2 in alignment:
			print (">"+s2.id, file=outhandle)
			if s2.id in trimdict:
				templistseq = list(s2.seq)
				for r in trimdict[s2.id]:
					if type(r) is int:
						templistseq[r] = "$"
						trimposcounter += 1
					else:
						templistseq[r[0]:r[1]] = "$"*endlen
						trimposcounter += endlen
				tempseq = "".join(templistseq)
			else:
				tempseq = s2.seq
			print (tempseq, file=outhandle)
	print (trimposcounter)
