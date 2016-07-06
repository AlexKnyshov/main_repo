from Bio import SeqIO
import glob
import os
import shutil
import csv
import sys
if len(sys.argv) == 6:
    blastfilearg = sys.argv[1]
    trif = sys.argv[2]
    ahefoldarg = sys.argv[3]
    evalue = float(sys.argv[4])
    opt = sys.argv[5]
else:
    print "FORMAT: python blparser.py [blastfile or folder] [asemblyfile or folder] [ahefolder] [evalue] [option: -n (normal), -s (extract only matched parts), -mn (multiple normal), -ms (multiple selected)]"
    print "EXAMPLE: python blparser.py blast.tab trinity.fas ./fasta 1e-40 -n"
    sys.exit()

output = {} #main dctionary
trans_d = {} #ahe loci and corresponding start and end of trans locus

#multi db option
if opt == "-mn" or "-ms":
    #prep
    if not os.path.exists ("./modified/"):
        os.makedirs("./modified") #creating folder if necessary
    else:
        shutil.rmtree("./modified/") #removing old files
        os.makedirs("./modified")

    #copy files
    print "copying files:"
    for x in glob.glob(ahefoldarg+"/*.fas"):
        locusfname = x.split("/")[-1]
        #print locusfname
        if not os.path.exists ("./modified/"+locusfname):
            prog = "copying "+str(locusfname)+"..."
            sys.stdout.write(prog+"\r")
            sys.stdout.flush()
            shutil.copy2(ahefoldarg+locusfname, "./modified")
    print ""
    #reading the blastfile
    blastlist = glob.glob(blastfilearg+"/*.blast")
    for b in blastlist:
        blastfile = open(b, "rU")
        reader = csv.reader(blastfile, delimiter='\t')
        currentkey = ""
        for row in reader:
            if currentkey != row[0]: ##new query
                if float(row[10]) <= evalue:
                    if row[0].split("//")[-1] in output: ##query is present
                        print "warning: the key exists", row[0].split("//")[-1]
                    else:
                        output[row[0].split("//")[-1]] = row[1] ##standart output
                        #output[AHE] = trans
                        #trans_d[AHE] = [start, finish] - for a given AHE, what was the start and end
                        if opt == "-s":
                            start = min(int(row[8]), int(row[9])) #secret opt
                            end = max(int(row[8]), int(row[9])) #secret opt
                            trans_d[row[0].split("//")[-1]] = [start, end] #secret opt
                        print row[0].split("//")[-1], row[1], row[2], row[10], row[11]
            else: ##same query
                if currentmatch != row[1] and float(row[10]) <= evalue:
                    print "warning: several matches detected", row[0].split("//")[-1].split(".tx_tm")[0], row[1], "delta is", currente-float(row[10])
            currentkey = row[0]
            currentmatch = row[1]
            currente = float(row[10])
        blastfile.close()
        count = int(len(output))
        print count, "targets found to be extracted"

        print "scanning the transcriptome..."
        #output2 = []
        warninglist = []
        inputf = SeqIO.parse(trif+"/"+b[:-6].split("/")[-1], "fasta")
        print "searching for contigs in:", b[:-6].split("/")[-1]
        c1 = 0
        extr_loci = []
        for seq in inputf:
            #print count
            if count == 0:
                print "search terminated"
                break
            elif seq.id in output.values(): #if contig is in ahe (was found as a blast hit)
                print "start"
                for x,y in output.items(): #checking all ahe (seach which ahe has it)
                    locusfname = x #AHE name
                    #y - trans locus name
                    #seq.id = trans locus name as well
                    if y == seq.id:
                        temp = seq.id
                        tempseq = seq.seq #secret opt
                        if opt == "-s":
                            seq.seq = seq.seq[trans_d[locusfname][0]:trans_d[locusfname][1]] #secret opt
                            print "s:", trans_d[locusfname][0], "e:", trans_d[locusfname][1], "trans end:", len(tempseq)
                            rhandle = open("./modified/"+locusfname, "r")
                            ali = SeqIO.parse(rhandle, "fasta")
                            for a in ali:
                                print "lenght AHE:", len(a), "length blast hit", len(seq.seq)
                                if float(len(seq.seq)) / len(a) < 0.8:
                                    print "Warning: blast hit is too short"
                                    warninglist.append(locusfname)
                                break
                            rhandle.close()
                        print "found", locusfname, y, "length:", len(seq.seq)
                        extr_loci.append(locusfname)
                        #print seq, len(seq) #debug
                        fhandle = open("./modified/"+locusfname, "a")
                        seq.id = b[:-6].split("/")[-1][:5]
                        seq.name =""
                        seq.description =""
                        if opt == "-ms" and locusfname not in warninglist:
                            SeqIO.write(seq, fhandle, "fasta")
                            c1 += 1
                        elif opt == "-mn":
                            SeqIO.write(seq, fhandle, "fasta")
                            c1 += 1
                        seq.id = temp
                        seq.seq = tempseq #secret opt
                        count -= 1
                        fhandle.close()
                    # else:
                    #     print x, y
        print c1, "loci extracted"
        count = int(len(output))
        print "files written: ", len(extr_loci)
        print "warning list:", len(warninglist)

else:
    #reading blast file
    print "reading blastfile...", blastfilearg
    blastfile = open(blastfilearg, "rU")
    reader = csv.reader(blastfile, delimiter='\t')
    currentkey = ""
    for row in reader:
        if currentkey != row[0]: ##new query
            if float(row[10]) <= evalue:
                if row[0].split("//")[-1] in output: ##query is present
                    print "warning: the key exists", row[0].split("//")[-1]
                else:
                    output[row[0].split("//")[-1]] = row[1] ##standart output
                    #output[AHE] = trans
                    #trans_d[AHE] = [start, finish] - for a given AHE, what was the start and end
                    if opt == "-s":
                        start = min(int(row[8]), int(row[9])) #secret opt
                        end = max(int(row[8]), int(row[9])) #secret opt
                        trans_d[row[0].split("//")[-1]] = [start, end] #secret opt
                    print row[0].split("//")[-1], row[1], row[2], row[10], row[11]
        else: ##same query
            if currentmatch != row[1] and float(row[10]) <= evalue:
                print "warning: several matches detected", row[0].split("//")[-1].split(".tx_tm")[0], row[1], "delta is", currente-float(row[10])
        currentkey = row[0]
        currentmatch = row[1]
        currente = float(row[10])
    blastfile.close()
    count = int(len(output))
    print count, "targets found to be extracted"

    #print output
    #print trans_d

    #scanning the transcriptomes
    if not os.path.exists ("./modified/"):
        os.makedirs("./modified") #creating folder if necessary
    else:
        shutil.rmtree("./modified/") #removing old files
        os.makedirs("./modified")

    #copy files
    print "copying files:"
    for x in glob.glob(ahefoldarg+"/*.fas"):
        locusfname = x.split("/")[-1]
        #print locusfname
        if not os.path.exists ("./modified/"+locusfname):
            prog = "copying "+str(locusfname)+"..."
            sys.stdout.write(prog+"\r")
            sys.stdout.flush()
            shutil.copy2(ahefoldarg+locusfname, "./modified")
    print ""

    print "scanning the transcriptome..."
    #output2 = []
    warninglist = []
    inputf = SeqIO.parse(trif, "fasta")
    print "searching for contigs in:", trif
    c1 = 0
    extr_loci = []
    for seq in inputf:
        #print count
        if count == 0:
            print "search terminated"
            break
        elif seq.id in output.values(): #if contig is in ahe (was found as a blast hit)
            print "start"
            for x,y in output.items(): #checking all ahe (seach which ahe has it)
                locusfname = x #AHE name
                #y - trans locus name
                #seq.id = trans locus name as well
                if y == seq.id:
                    temp = seq.id
                    tempseq = seq.seq #secret opt
                    if opt == "-s":
                        seq.seq = seq.seq[trans_d[locusfname][0]:trans_d[locusfname][1]] #secret opt
                        print "s:", trans_d[locusfname][0], "e:", trans_d[locusfname][1], "trans end:", len(tempseq)
                        rhandle = open("./modified/"+locusfname, "r")
                        ali = SeqIO.parse(rhandle, "fasta")
                        for a in ali:
                            print "lenght AHE:", len(a), "length blast hit", len(seq.seq)
                            if float(len(seq.seq)) / len(a) < 0.8:
                                print "Warning: blast hit is too short"
                                warninglist.append(locusfname)
                            break
                        rhandle.close()
                    print "found", locusfname, y, "length:", len(seq.seq)
                    extr_loci.append(locusfname)
                    #print seq, len(seq) #debug
                    fhandle = open("./modified/"+locusfname, "a")
                    seq.id = trif.split("/")[-1][:5]
                    seq.name =""
                    seq.description =""
                    if opt == "-s" and locusfname not in warninglist:
                        SeqIO.write(seq, fhandle, "fasta")
                        c1 += 1
                    elif opt == "-n":
                        SeqIO.write(seq, fhandle, "fasta")
                        c1 += 1
                    seq.id = temp
                    seq.seq = tempseq #secret opt
                    count -= 1
                    fhandle.close()
                # else:
                #     print x, y
    print c1, "loci extracted"
    count = int(len(output))
    print "files written: ", len(extr_loci)
    print "warning list:", len(warninglist)
print "done"