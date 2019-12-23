import os
import sys
from Bio import AlignIO

infile = open(sys.argv[1], "r") #input alignment
window = int(sys.argv[2]) #detection window -- tested 20
step = int(sys.argv[3]) #stepsize -- tested 20
thresh = int(sys.argv[4]) #how many bases in window shoudl be bad - probs need half bad? -- tested 10
misdata = float(sys.argv[5]) #how many bases in window shoudl be bad - probs need half bad? -- tested 10
misdata_seq = float(sys.argv[6]) #how many bases in window shoudl be bad - probs need half bad? -- tested 10
trimsm = True

inputalignment = AlignIO.read(infile, "fasta")
infile.close()

outname = sys.argv[1].split("/")[-1]

length = inputalignment.get_alignment_length()

seq_order = []
for seq in inputalignment:
	seq_order.append(seq.id)

# print len(seq_order)

outdict = {}
outdict2 = {}
outdict3 = {}
outdict4 = {}
outdict5 = {}
for tx in inputalignment:
	outdict[tx.id] = [] #basemask
	outdict2[tx.id] = []
	outdict3[tx.id] = []
	outdict4[tx.id] = []
	outdict5[tx.id] = []
	#print outdict[tx.id][0]

#build profile at 40% threshold
# nucl_profile = ""
gap_profile = ""
for x in range(length):
	col = inputalignment[:,x].lower()
	a = col.count('a')
	t = col.count('t')
	g = col.count('g')
	c = col.count('c')
	tot = a+t+g+c
	if tot > 0:
		ok_bases = set()
		ok_bases.add("-")
		if float(a) / tot > 0.4:
			ok_bases.add("a")
		if float(t) / tot > 0.4:
			ok_bases.add("t")
		if float(g) / tot > 0.4:
			ok_bases.add("g")
		if float(c) / tot > 0.4:
			ok_bases.add("c")
		#print ok_bases
		if float(col.count("-")) / (tot+col.count("-")) > misdata:
			gap_profile += "-"
		else:
			gap_profile += "n"
		for seq in range(len(col)):
			#print inputalignment[seq].id
			if col[seq] not in ok_bases:
				outdict[inputalignment[seq].id].append("n")
			else:
				outdict[inputalignment[seq].id].append(col[seq])
			outdict2[inputalignment[seq].id].append(col[seq])
			outdict3[inputalignment[seq].id].append(col[seq])
			outdict4[inputalignment[seq].id].append(col[seq])
			outdict5[inputalignment[seq].id].append(col[seq])
	else:
		gap_profile += "-"
		for seq in range(len(col)):
			outdict[inputalignment[seq].id].append("-")
		# outdict[inputalignment[seq].id].append(col[seq])
		# outdict2[inputalignment[seq].id].append(col[seq])
			outdict2[inputalignment[seq].id].append(col[seq])
			outdict3[inputalignment[seq].id].append(col[seq])
			outdict4[inputalignment[seq].id].append(col[seq])
			outdict5[inputalignment[seq].id].append(col[seq])
	# print len(gap_profile), x+1
	# if len(gap_profile) != x+1:
	# 	sys.exit()

# outf = open("basemask"+outname, "w")
# mapped = {v: i for i, v in enumerate(seq_order)}

# for item, val in sorted(outdict.items(), key=lambda x: mapped[x[0]]):
# #for item, val in outdict.items():
# 	print >> outf, ">"+item
# 	print >> outf, "".join(val)
# outf.close()

#per seq masking
gap_profile_2 = [0 for x in range(length)]
for seqid, seqval in outdict.items():
	#bad bases masking
	i = 0
	while i < length:
		if i+window<=length:
			end = i+window
		else:
			end = length
		# string = "".join(seqval[i:end])
		# gap_prop = float(gap_profile[i:end].count("-")) / len(gap_profile[i:end])
		gap_index = []
		nucl_index = []
		bad_num = 0
		#print i, end, length, len(gap_profile)
		for pos in range(i, end):
			if gap_profile[pos] == "-": #gap position
				gap_index.append(pos)
			elif gap_profile[pos] == "n": #nucl pos
				if seqval[pos] != "-":
					nucl_index.append(pos)
				if seqval[pos] == "n":# or seqval[pos] == "-": #count bad bases
					bad_num += 1
		#gap positions trimming
		if len(gap_index) > 0:
			for pos in gap_index:
		 		if seqval[pos] != "-":
		 			outdict2[seqid][pos] = "n"
		#bad positions trimming
		if len(nucl_index) > 0.25 * window: #0: #trim only if there are enough bases, like a quarter of window at least
			if float(bad_num) * window / len(nucl_index) > thresh: #if total is bad, mask; rescaled to nucl amount
				for pos in nucl_index:
			 		if seqval[pos] != "-":
			 			outdict2[seqid][pos] = "n"
			# else: #unmask
			# 	for pos in nucl_index:
			#  		outdict2[seqid][pos] = outdict4[seqid][pos]
			
		i += step

	#short seq masking
	i = 0
	while i < length:
		if i+window<=length:
			end = i+window
		else:
			end = length
		# print int(window*0.25)
		if outdict2[seqid][i:end].count("n") > thresh:
			for pos in range(i, end):
		 		if outdict2[seqid][pos] != "-":
		 			outdict3[seqid][pos] = "n"
					
		#elif outdict2[seqid][i:(i+2)].count("-") >= 2 and outdict2[seqid][(end-2):end].count("-") >= 2:
		elif trimsm and outdict2[seqid][i:(i+2)].count("-") == 2 and outdict2[seqid][(end-2):end].count("-") == 2:
			for pos in range(i, end):
		 		if outdict2[seqid][pos] != "-":
		 			outdict3[seqid][pos] = "n"

		elif i == 0 and outdict2[seqid][i:end].count("n")+outdict2[seqid][i:end].count("-") > thresh:
			for pos in range(i, end):
		 		if outdict2[seqid][pos] != "-":
		 			outdict3[seqid][pos] = "n"

		elif end == length and outdict2[seqid][i:end].count("n")+outdict2[seqid][i:end].count("-") > thresh:
			for pos in range(i, end):
		 		if outdict2[seqid][pos] != "-":
		 			outdict3[seqid][pos] = "n"

		# else: #unmask
		# 	for pos in nucl_index:
		# 		outdict3[seqid][pos] = outdict4[seqid][pos]
		i += step

	for pos in range(length):
		if outdict3[seqid][pos] == "n" or outdict3[seqid][pos] == "-":
			gap_profile_2[pos] += 1
		else:
			gap_profile_2[pos] += 0
	# # 	#back pass
	# i = length-1
	# while i >= step-1:
	# 	#print i, step, i-window
	# 	if i-window>=0:
	# 		start = i-window
	# 	else:
	# 		start = 0
	# 	if outdict2[seqid][i:end].count("n")+outdict2[seqid][i:end].count("-") > thresh:
	# 		for pos in range(i, end):
	# 	 		if outdict2[seqid][pos] != "-":
	# 	 			outdict5[seqid][pos] = "n"
	# 	else: #unmask
	# 		for pos in nucl_index:
	# 			outdict5[seqid][pos] = outdict4[seqid][pos]
	# 	i -= step




outf2 = open("masked"+outname, "w")
outf3 = open("trimmed_"+outname, "w")

for item1, val1 in sorted(outdict3.items(), key=lambda x: mapped[x[0]]):
	print >> outf2, ">"+item1
	print >> outf2, "".join(val1)
	# print >> outf2, ">"+item1+"out2"
	# print >> outf2, "".join(outdict2[item1])
	# print >> outf2, ">"+item1+"out3"
	# print >> outf2, "".join(outdict3[item1])


	new_val = ""
	for pos in range(length):
		if float(gap_profile_2[pos]) / len(outdict3) < misdata:
			if val1[pos] == "n":
				new_val += "-"
			else:
				new_val += val1[pos]
	if float(new_val.count("-")) / len(new_val) < misdata_seq:
		print >> outf3, ">"+item1
		print >> outf3, new_val
# print >> outf2, ">profile"
# print >> outf2, nucl_profile
print >> outf2, ">gap_profile"
print >> outf2, gap_profile
outf2.close()

outf3.close()

###version 2
# import os
# import sys
# from Bio import AlignIO

# infile = open(sys.argv[1], "r") #input alignment
# window = int(sys.argv[2]) #detection window -- tested 20
# step = int(sys.argv[3]) #stepsize -- tested 20
# thresh = int(sys.argv[4]) #how many bases in window shoudl be bad - probs need half bad? -- tested 10
# trimsm = True

# inputalignment = AlignIO.read(infile, "fasta")
# infile.close()

# outname = sys.argv[1].split("/")[-1]

# length = inputalignment.get_alignment_length()
# seq_order = []
# for seq in inputalignment:
# 	seq_order.append(seq.id)

# outdict = {}
# outdict2 = {}
# outdict3 = {}
# outdict4 = {}
# for tx in inputalignment:
# 	outdict[tx.id] = []
# 	outdict2[tx.id] = []
# 	outdict3[tx.id] = []
# 	outdict4[tx.id] = []
# 	#print outdict[tx.id][0]

# #build profile at 40% threshold
# # nucl_profile = ""
# gap_profile = ""
# for x in range(length):
# 	col = inputalignment[:,x].lower()
# 	a = col.count('a')
# 	t = col.count('t')
# 	g = col.count('g')
# 	c = col.count('c')
# 	tot = a+t+g+c
# 	if tot > 0:
# 		ok_bases = set()
# 		ok_bases.add("-")
# 		if float(a) / tot > 0.4:
# 			ok_bases.add("a")
# 		if float(t) / tot > 0.4:
# 			ok_bases.add("t")
# 		if float(g) / tot > 0.4:
# 			ok_bases.add("g")
# 		if float(c) / tot > 0.4:
# 			ok_bases.add("c")
# 		#print ok_bases
# 		if float(col.count("-")) / (tot+col.count("-")) > 0.7:
# 			gap_profile += "-"
# 		else:
# 			gap_profile += "n"
# 		for seq in range(len(col)):
# 			#print inputalignment[seq].id
# 			if col[seq] not in ok_bases:
# 				outdict[inputalignment[seq].id].append("n")
# 			else:
# 				outdict[inputalignment[seq].id].append(col[seq])
# 			outdict2[inputalignment[seq].id].append(col[seq])
# 			outdict3[inputalignment[seq].id].append(col[seq])
# 			outdict4[inputalignment[seq].id].append(col[seq])
# 	else:
# 		# outdict[inputalignment[seq].id].append(col[seq])
# 		# outdict2[inputalignment[seq].id].append(col[seq])
# 		outdict[inputalignment[seq].id].append("n")
# 		outdict2[inputalignment[seq].id].append(col[seq])
# 		outdict3[inputalignment[seq].id].append(col[seq])
# 		outdict4[inputalignment[seq].id].append(col[seq])

# outf = open("basemask"+outname, "w")
# mapped = {v: i for i, v in enumerate(seq_order)}

# for item, val in sorted(outdict.items(), key=lambda x: mapped[x[0]]):
# #for item, val in outdict.items():
# 	print >> outf, ">"+item
# 	print >> outf, "".join(val)
# outf.close()

# #per seq masking
# for seqid, seqval in outdict.items():
# 	i = 0
# 	while i < length:
# 		if i+window<=length:
# 			end = i+window
# 		else:
# 			end = length
# 		# string = "".join(seqval[i:end])
# 		# gap_prop = float(gap_profile[i:end].count("-")) / len(gap_profile[i:end])
# 		gap_index = []
# 		nucl_index = []
# 		bad_num = 0
# 		for pos in range(i, end):
# 			if gap_profile[pos] == "-": #gap position
# 				gap_index.append(pos)
# 			elif gap_profile[pos] == "n": #nucl pos
# 				nucl_index.append(pos)
# 				if seqval[pos] == "n" or seqval[pos] == "-": #count bad bases
# 					bad_num += 1
# 		if len(gap_index) > 0:
# 			for pos in gap_index:
# 		 		if seqval[pos] != "-":
# 		 			outdict2[seqid][pos] = "n"
# 		if len(nucl_index) > 0:
# 			if float(bad_num) * window / len(nucl_index) > thresh: #if total is bad, mask; rescaled to nucl amount
# 				for pos in nucl_index:
# 			 		if seqval[pos] != "-":
# 			 			outdict2[seqid][pos] = "n"
# 			else: #unmask
# 				for pos in nucl_index:
# 			 		outdict2[seqid][pos] = outdict4[seqid][pos]
			
# 		# if gap_prop < 0.1 and string.count("n")+string.count("-") > thresh: #internal zone - check for orphan
# 		# 	for pos in range(i,end):
# 		# 		if seqval[pos] != "-":
# 		# 			outdict2[seqid][pos] = "n"
# 		# elif gap_prop > 0.9 and string.count("n") > thresh:
# 		# 	for pos in range(i,end):
# 		# 		outdict2[seqid][pos] = "n"
# 		# 		#+string.count("a")+string.count("t")+string.count("c")+string.count("g")
# 		# elif gap_prop < 0.1 and string.count("n")+string.count("-") <= thresh:
# 		# 	for pos in range(i,end):
# 		# 		outdict2[seqid][pos] = outdict4[seqid][pos]
		
# 		i += step
# 	#back pass
# 	i = length-1
# 	while i >= step:
# 		#print i, step, i-window
# 		if i-window>=0:
# 			start = i-window
# 		else:
# 			start = 0
# 		# string = "".join(seqval[start:i])
# 		# gap_prop = float(gap_profile[start:i].count("-")) / len(gap_profile[start:i])
# 		# if gap_prop < 0.1 and string.count("n")+string.count("-") > thresh: #internal zone - check for orphan
# 		# 	for pos in range(start,i):
# 		# 		if seqval[pos] != "-":
# 		# 		 	outdict3[seqid][pos] = "n"
# 		# elif gap_prop > 0.9 and string.count("n") > thresh:
# 		# 	for pos in range(start,i):
# 		# 		outdict3[seqid][pos] = "n"
# 		# elif gap_prop < 0.1 and string.count("n")+string.count("-") <= thresh:
# 		# 	for pos in range(i,end):
# 		# 		outdict3[seqid][pos] = outdict4[seqid][pos]

# 		gap_index = []
# 		nucl_index = []
# 		bad_num = 0
# 		for pos in range(start,i):
# 			if gap_profile[pos] == "-": #gap position
# 				gap_index.append(pos)
# 			elif gap_profile[pos] == "n": #nucl pos
# 				nucl_index.append(pos)
# 				if seqval[pos] == "n" or seqval[pos] == "-": #count bad bases
# 					bad_num += 1
# 		if len(gap_index) > 0:
# 			for pos in gap_index:
# 		 		if seqval[pos] != "-":
# 		 			outdict3[seqid][pos] = "n"
# 		if len(nucl_index) > 0:
# 			if float(bad_num) * window / len(nucl_index) > thresh: #if total is bad, mask
# 				for pos in nucl_index:
# 			 		if seqval[pos] != "-":
# 			 			outdict3[seqid][pos] = "n"
# 			else: #unmask
# 				for pos in nucl_index:
# 			 		outdict3[seqid][pos] = outdict4[seqid][pos]
			
# 		i -= step
# 	#final pass
# 	for p in range(length):
# 		if outdict3[seqid][p] == outdict2[seqid][p] == "n":
# 			outdict4[seqid][p] = "n"

# 		#print outdict2[seqid][pos], outdict3[seqid][pos]

# 				#+string.count("a")+string.count("t")+string.count("c")+string.count("g")
# 		# elif gap_prop < 0.1 and string.count("n") <= thresh:
# 		# 	for pos in range(i,end):
# 		# 		outdict2[seqid][pos] = outdict3[seqid][pos]
# 		# elif gap_prop <= 0.9 and gap_prop >= 0.1 and string.count("n") > thresh:
# 		# 	for pos in range(i,end):
# 		# 		outdict2[seqid][pos] = "n"

# outf2 = open("masked"+outname, "w")
# #mapped = {v: i for i, v in enumerate(seq_order)}

# for item1, val1 in sorted(outdict4.items(), key=lambda x: mapped[x[0]]):
# 	print >> outf2, ">"+item1
# 	print >> outf2, "".join(val1)
# 	# print >> outf2, ">"+item1+"out2"
# 	# print >> outf2, "".join(outdict2[item1])
# 	# print >> outf2, ">"+item1+"out3"
# 	# print >> outf2, "".join(outdict3[item1])
# # print >> outf2, ">profile"
# # print >> outf2, nucl_profile
# print >> outf2, ">gap_profile"
# print >> outf2, gap_profile
# outf2.close()
