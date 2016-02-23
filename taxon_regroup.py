import sys
import csv
from Bio import SeqIO
from Bio import AlignIO


csvname = sys.argv[1]
aligname = sys.argv[2]

db = {}

Enicos = []
Ceratos = []
Dipsos = []
Hypsos = []
Ogeriinae = []
Schizos = []
Others = []

Hypsolist = ["Duonota", "Glyptocombus", "Hypselosoma", "Lativena", "Macromannus", "Pateena", "Rectilamina", "Williamsocoris", "Ommatides", "Hypselosomatinae", "NrHypselosoma"]
Ogerolist = ["Chinannus", "Kaimon", "Kokeshia", "Ogeria", "nr Kokeshia", "nr Ogeria", "nr Pachyplagioides"]
Schizolist = ["Ceratocomboides", "Corixidea", "Dundonannus", "Hoplonannus", "Membracioides", "Nannocoris", "Pinochius", "Ptenidiophyes", "Schizoptera", "Voccoroda", "Voragocoris", "nr Ceratocomboides", "nr Corixidea", "nr Dundonannus", "nr Hoplonannus", "nr Membracioides", "nr Ptenidiophyes", "nr Semanganannus", "nr Semangananus", "Nannodictyus"]

with open(csvname, "rb") as csvfile:
	reader = csv.reader(csvfile)
	taxonset = set()
	for row in reader:
		ED = int(row[0].split("_")[1])
		infra = row[2]
		family = row[3]
		genus = row[5]
		#print ED, infra	,family, genus
		if infra == "Enicocephalomorpha":
			db[(str(ED))] = "Enico"
		elif family == "Ceratocombidae":
			db[(str(ED))] = "Cerato"
		elif family == "Dipsocoridae":
			db[(str(ED))] = "Dipso"
		elif family == "Schizopteridae":
			if genus in Hypsolist:
				db[(str(ED))] = "Hypso"
			elif genus in Ogerolist:
				db[(str(ED))] = "Ogeriinae"
			elif genus in Schizolist:
				db[(str(ED))] = "Schizo"
			else:
				db[(str(ED))] = "Other"
				#taxonset.add(genus)
		else:
			print "do nothing"
# for x in sorted(taxonset):
# 	print x
input_handle = open(aligname, "rU")
alignments = SeqIO.parse(input_handle, "fasta")
index = 0
for seq in alignments:
	EDfasta = seq.id.split("_")[-1]
	genusfasta = seq.id.split("_")[0]
	if EDfasta in db:
		if db[EDfasta] == "Enico":
			Enicos.append(index)

		if db[EDfasta] == "Cerato":
			Ceratos.append(index)

		if db[EDfasta] == "Dipso":
			Dipsos.append(index)

		if db[EDfasta] == "Hypso":
			Hypsos.append(index)

		if db[EDfasta] == "Ogeriinae":
			Ogeriinae.append(index)

		if db[EDfasta] == "Schizo":
			Schizos.append(index)

		if db[EDfasta] == "Other":
			if genusfasta in Hypsolist:
				Hypsos.append(index)
			elif genusfasta in Ogerolist:
				Ogeriinae.append(index)
			elif genusfasta in Schizolist:
				Schizos.append(index)
			else:
				Others.append(index)
	else:
		print seq.id
		Others.append(index)
	index += 1
taxa = []
taxa.append(Enicos)
taxa.append(Ceratos)
taxa.append(Dipsos)
taxa.append(Hypsos)
taxa.append(Ogeriinae)
taxa.append(Schizos)
taxa.append(Others)
#Dipsos = []
#Hypsos = []
#Ogeriinae = []
#Schizos = []
#Others = []
#print taxa
input_handle.seek(0)
alignments = AlignIO.read(input_handle, "fasta")
#print alignments
outfile = open(aligname+".txt", "w")
group = 1
for g in taxa:
	Goutfile = open(aligname+str(group)+".fas", "w")
	for t in g:
		print >> outfile, ">"+alignments[t].id
		print >> outfile, alignments[t].seq
		print >> Goutfile, ">"+alignments[t].id
		print >> Goutfile, alignments[t].seq 
	Goutfile.close()
	group += 1
	# if seq.id in d:
	# 	d[seq.id].append(infile.split("/")[-1])
	# else:
	# 	d[seq.id] = []
	# 	d[seq.id].append(infile.split("/")[-1])
outfile.close()