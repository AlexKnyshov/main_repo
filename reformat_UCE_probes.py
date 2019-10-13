from Bio import SeqIO
import sys


if len(sys.argv) >= 3:
	infname = sys.argv[1]
	taxname = sys.argv[2]

else:
	print "FORMAT: argument1 = probe file in fasta format, argument2 = taxon name"
	print "EXAMPLE: ./fasta cimlec1"
	print "output is written to [taxon name]_probe_regions.fas"
	sys.exit()

newseqs = {} #probes_locus = [start, end, seq]
with open(infname) as infhandle:
	seqs = SeqIO.parse(infhandle, "fasta")
	for seq in seqs:
		seqlist = seq.description.split("|")[-1].split(",")
		seqdict = {}
		for f in seqlist:
			fsplit = f.split(":")
			seqdict[fsplit[0]] = fsplit[1]
		if seqdict['probes-source'] == taxname:
			if seqdict['probes-locus'] in newseqs.keys():
				if newseqs[seqdict['probes-locus']][0] < int(seqdict['probes-local-start']):
					lendif = int(seqdict['probes-local-end']) - newseqs[seqdict['probes-locus']][1]
					newseqs[seqdict['probes-locus']][2] = newseqs[seqdict['probes-locus']][2]+str(seq.seq[len(seq.seq)-lendif:])
					newseqs[seqdict['probes-locus']][1] = int(seqdict['probes-local-end'])
				else:
					lendif = newseqs[seqdict['probes-locus']][0] - int(seqdict['probes-local-start'])
					newseqs[seqdict['probes-locus']][2] = str(seq.seq[:lendif])+newseqs[seqdict['probes-locus']][2]
					newseqs[seqdict['probes-locus']][0] = int(seqdict['probes-local-start'])
			else:
				newseqs[seqdict['probes-locus']] = [int(seqdict['probes-local-start']),int(seqdict['probes-local-end']),str(seq.seq)]
with open(taxname+"_probe_regions.fas", "w") as outh:
	for key, val in newseqs.items():
		print >> outh, ">"+key
		print >> outh, val[2]