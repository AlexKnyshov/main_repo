from Bio import AlignIO
import sys
import glob
import os
import shutil

inputfolder = sys.argv[1]
if len(sys.argv)>2:
	exclusion_file = sys.argv[2]

exclusion_list = []
if len(sys.argv) == 3:
	print "reading exclusion list..."
	exfile = open(exclusion_file, "r")
	for line in exfile:
		l = line.strip()
		exclusion_list.append(l)
	exfile.close()
	print "read", len(exclusion_list), "records"


files = glob.glob(inputfolder+"/*.fas")
count = 0
counter = 0
misdata = 0
totaldata = 0
locilist = []
for f in files:
	count +=1
	if f.split("/")[-1] not in exclusion_list:
		counter += 1
		locilist.append(f)
	#for seq_record in SeqIO.parse(f, "fasta"):
	#	counter +=1

	# ali = AlignIO.read(f, "fasta")
	# print f, len(ali)
	# print f, len(ali)#counter
	# if len(ali) > 8:
	# 	counter += 1
	# 	print f, "good", len(ali)
	# 	locilist.append(f)

	# for seq in ali:
	# 	if seq.id == "Cryptostemma_sp_Peru_249":
	# 		counter += 1
	# 		print f, "good", len(ali)
	# 		locilist.append(f)

print "Good %i out of total %i records" % (counter, count)

#copy files
if not os.path.exists ("./subset2/"):
    os.makedirs("./subset2") #creating folder if necessary
else:
    shutil.rmtree("./subset2/") #removing old files
    os.makedirs("./subset2")

print "copying files:"
for x in locilist:
    locusfname = x.split("/")[-1]
    #print locusfname
    if not os.path.exists ("./subset2/"+locusfname):
        prog = "copying "+str(locusfname)+"..."
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        shutil.copy2(inputfolder+locusfname, "./subset2")
print "done"