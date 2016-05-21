from Bio import SeqIO
import sys
import glob
import os
import shutil

f = sys.argv[1]
ext = sys.argv[2]
files = glob.glob(f+"/*"+ext)


# d = {}
# trimlist = {}
# trimlist["I13432_ED_4993_Hemiptera_Dipsocoridae_Cryptostemma_sp_seq1"] = "Cryptostemma_sp_Peru_249"
# trimlist["I13433_ED_2045_Hemiptera_Ceratocombidae_Kvamula_sp_seq1"] = "cf_Kvamula_or_Seychellesanus_sp_Madagascar_2043"
# trimlist["I13434_ED_2660_Hemiptera_Schizopteridae_Williamsocoris_sp_seq1"] = "Williamsocoris_sp_Trinidad_2660"
# trimlist["I13435_ED_4258_Hemiptera_Schizopteridae_Nannocoris_sp_seq1"] = "Nannocoris_sp_4258"
# trimlist["I13436_ED_1692_Hemiptera_Schizopteridae_Kokeshia_sp_seq1"] = "Kokeshia_sp_Thailand_1409"
# trimlist["I13437_ED_4257_Hemiptera_Schizopteridae_Chinannus_sp_seq1"] = "Chinannus_monteverdensis_4257"
# trimlist["I13438_ED_2192_Hemiptera_Schizopteridae_Dundonannus_sp_seq1"] = "Dundonannus_sp_2190"
# trimlist["I13439_ED_6303_Hemiptera_Schizopteridae_Schizoptera_sp_seq1"] = "Schizoptera_sp_6303"

# #trimlist["I13445_RCW_4101_Hemiptera_Phymatinae_Phymata_pacifica_seq1"] = "Phymata"

# print "input"
# for infile in files:
# 	input_handle = open(infile, "rU")
# 	alignments = SeqIO.parse(input_handle, "fasta")
# 	print "read", infile
# 	outhandle = open(infile+".tx_tm.fas", "w")
# 	for seq in alignments:
# 		#if seq.id in trimlist:
# 		if seq.id in trimlist.values():
# 			print >> outhandle, ">"+trimlist[seq.id]
# 			print >> outhandle, seq.seq
# print "done"

translist = ["EnspE", "Lican", "CespC", "Beflu", "GAYI0"]
locilist = set()
print "input"
for infile in files:
	input_handle = open(infile, "rU")
	alignments = SeqIO.parse(input_handle, "fasta")
	print "read", infile
	#outhandle = open(infile+".tx_tm.fas", "w")
	for seq in alignments:
		#if seq.id in trimlist:
		if seq.id in translist:
			locilist.add(infile)
print locilist
print len(locilist)

#copy files
if not os.path.exists ("./subset/"):
    os.makedirs("./subset") #creating folder if necessary
else:
    shutil.rmtree("./subset/") #removing old files
    os.makedirs("./subset")

print "copying files:"
for x in locilist:
    locusfname = x.split("/")[-1]
    #print locusfname
    if not os.path.exists ("./subset/"+locusfname):
        prog = "copying "+str(locusfname)+"..."
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        shutil.copy2(f+locusfname, "./subset")
print "done"
