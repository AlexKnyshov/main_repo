from Bio import SeqIO
import glob
import os
import shutil
import csv
import sys
import math
if len(sys.argv) == 6:
    blastfilearg = sys.argv[1]
    trif = sys.argv[2]
    ahefoldarg = sys.argv[3]
    evalue = float(sys.argv[4])
    opt = sys.argv[5]
else:
    print "FORMAT: python alexblparser.py [blastfile or folder] [asemblyfile or folder] [ahefolder] [evalue] [option: -n (normal), -s (extract only matched parts), -ss (short are discarded), -mn (multiple normal), -ms (multiple selected), -mss(short are discarded), -me (extended ms option)]"
    print "EXAMPLE: python alexblparser.py blast.tab trinity.fas ./fasta 1e-40 -n"
    sys.exit()

print "option", opt, "selected"
#multi db option
if opt == "-mn" or opt == "-ms" or opt == "-mss" or opt == "-me":
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
    #debug vars
    number = 0
    numberset = set()
    totalloci = 0
    #parsing blast files
    for b in blastlist:
        print "processing", b
        output = {} #main dctionary
        trans_d = {} #ahe loci and corresponding start and end of trans locus
        e_opt = {} #for a ahe (assuming the first hit is extracted), store trans and ahe data
        #e_opt[AHE] = [transf, transr, transb, ahef, aher, aheb]
        blastfile = open(b, "rU")
        reader = csv.reader(blastfile, delimiter='\t')
        currentkey = ""
        best_eval = 0.0
        best_hit = 0
        for row in reader:
            if currentkey != row[0]: ##new query
                if float(row[10]) <= evalue:
                    if row[0].split("//")[-1] in output: ##query is present
                        print "warning: the key exists", row[0].split("//")[-1]
                    else:
                        output[row[0].split("//")[-1]] = row[1] ##standart output
                        best_eval = float(row[10])
                        best_hit = float(row[11])
                        #output[AHE] = trans
                        #trans_d[AHE] = [start, finish] - for a given AHE, what was the start and end
                        if opt == "-ms" or opt == "-mss" or opt == "-me":
                            start = min(int(row[8]), int(row[9])) #secret opt
                            end = max(int(row[8]), int(row[9])) #secret opt
                            trans_d[row[0].split("//")[-1]] = [start, end] #secret opt
                            #check trans
                            if int(row[8]) < int(row[9]):
                                print int(row[8]), int(row[9]), "forward"
                                transb = True
                                #transf = int(row[8])
                                #transr = int(row[9])
                                #print "AHE f:", int(row[6]), int(row[9]), "forward"
                            else:
                                print int(row[8]), int(row[9]), "reverse"
                                transb = False
                            transf = int(row[8])
                            transr = int(row[9])
                            #check query
                            if int(row[6]) < int(row[7]):
                                print int(row[6]), int(row[7]), "forward"
                                aheb = True
                                #transf = int(row[8])
                                #transr = int(row[9])
                                #print "AHE f:", int(row[6]), int(row[9]), "forward"
                            else:
                                print int(row[7]), int(row[6]), "reverse"
                                aheb = False
                            ahef = int(row[6])
                            aher = int(row[7])
                            #store the data:
                            e_opt[row[0].split("//")[-1]] = [transf, transr, transb, ahef, aher, aheb]
                        print "--------------------------------------------------------"
                        print row[0].split("//")[-1], row[1], row[2], row[10], row[11]
            else: ##same query
                if currentmatch != row[1] and float(row[10]) <= evalue:
                    if row[0].split("//")[-1] in output: 
                        #adopted for trinity assemblies
                        if row[1].split("_c")[0] == output[row[0].split("/")[-1]].split("_c")[0]:
                            print "warning: several isoforms detected", row[1], row[10], row[11]
                        else:
                            print "warning: several matches detected", row[1], row[10], row[11], "hit is", round(float(row[11])/best_hit*100), "%"#"; ratio with best_eval is"#, round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                            if round(float(row[11])/best_hit*100) > 90:#round(math.log(float(row[10]))/math.log(best_eval)*100) > 90 or :
                                number += 1
                                numberset.add(row[0])
                    else:
                        output[row[0].split("/")[-1]] = row[1] ##standart output
                        print "ADDITION"
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
                        if opt == "-ms" or opt == "-mss" or opt == "-me":
                            rhandle = open("./modified/"+locusfname, "r")
                            ali = SeqIO.parse(rhandle, "fasta")
                            recs = list(ali)
                            #e option:
                            #e_opt[AHE] = [transf, transr, transb, ahef, aher, aheb]
                            #forward ahe:
                            #print "E opt"
                            if e_opt[locusfname][5]:#aheb: #ahe is forward
                                #print "AHE forw"
                                gap = e_opt[locusfname][3]#ahef #ahe start gap
                                if e_opt[locusfname][2]:#transb: #trans is forward
                                    ##need to check if trans has start - OK
                                    if e_opt[locusfname][0] - gap >= 0:
                                        gapt = e_opt[locusfname][0] - gap#transf - gap #this is the corrected trans start
                                    else:
                                        gapt = 0 # otherwise start from start
                                else:
                                    ##need to check if trans is long enough - OK
                                    if e_opt[locusfname][0] + gap <= len(tempseq):
                                        gapt = e_opt[locusfname][0] + gap#transf + gap #this is the corrected trans start
                                    else:
                                        gapt = len(tempseq) # otherwise start from the end
                                #check the ahe end:
                                gaprev = len(recs[0].seq) - e_opt[locusfname][4]#ali.get_alignment_length() - e_opt[locusfname][0][4]
                                #if gaprev >= 0:#aher #some AHE left
                                if e_opt[locusfname][2]:#transb: #trans is forward
                                    #check the trans end
                                    if e_opt[locusfname][1] + gaprev <= len(tempseq): #check that trans is long enough
                                        gaptrev = e_opt[locusfname][1] + gaprev#transf - gap #this is the corrected trans start
                                    else:
                                        gaptrev = len(tempseq)
                                else:
                                    ##need to check if trans has start
                                    if e_opt[locusfname][1] - gaprev >= 0:
                                        gaptrev = e_opt[locusfname][1] - gaprev#transf + gap #this is the corrected trans start
                                    else:
                                        gaptrev = 0
                            #reverse AHE:
                            else: #ahe is reverse
                                #print "AHE rev"
                                gap = len(recs[0].seq) - e_opt[locusfname][3]#ali.get_alignment_length() - e_opt[locusfname][0][3]#ahef #ahe start gap
                                if e_opt[locusfname][2]:#transb: #trans is forward
                                    ##need to check if trans has start - OK
                                    if e_opt[locusfname][0] - gap >= 0:
                                        gapt = e_opt[locusfname][0] - gap#transf - gap #this is the corrected trans start
                                    else:
                                        gapt = 0
                                else:
                                    ##need to check if trans is long enough - OK
                                    if e_opt[locusfname][0] + gap <= len(tempseq):
                                        gapt = e_opt[locusfname][0] + gap#transf + gap #this is the corrected trans start
                                    else:
                                        gapt = len(tempseq)
                                #check the ahe end:
                                gaprev = e_opt[locusfname][4] # ahe end gap
                                #if gaprev >= 0:#aher #some AHE left
                                if e_opt[locusfname][2]:#transb: #trans is forward
                                    #check the trans end
                                    if e_opt[locusfname][1] + gaprev <= len(tempseq): #check that trans is long enough
                                        gaptrev = e_opt[locusfname][1] + gaprev#transf - gap #this is the corrected trans start
                                    else:
                                        gaptrev = len(tempseq)
                                else:
                                    ##need to check if trans has start
                                    if e_opt[locusfname][1] - gaprev >= 0:
                                        gaptrev = e_opt[locusfname][1] - gaprev#transf + gap #this is the corrected trans start
                                    else:
                                        gaptrev = 0
                            print "BLAST:", "transf:",e_opt[locusfname][0], "transr:",e_opt[locusfname][1], "transb:",e_opt[locusfname][2], "ahef:",e_opt[locusfname][3], "aher:",e_opt[locusfname][4], "aheb:",e_opt[locusfname][5]
                            print "LOCUS", "ahe:", len(recs[0].seq), "trans:", len(seq.seq)
                            if opt == "-me":
                                #TEST
                                if e_opt[locusfname][2]:
                                    print "STATS:", e_opt[locusfname][2], "start:", gapt,"stop",  gaptrev
                                    seq.seq = seq.seq[gapt:gaptrev]
                                else:
                                    print "STATS:", e_opt[locusfname][2], "start:", gaptrev,"stop",  gapt
                                    seq.seq = seq.seq[gaptrev:gapt]
                                    seq = seq.reverse_complement()
                            elif opt == "-ms" or opt == "-mss":
                                #old option
                                seq.seq = seq.seq[trans_d[locusfname][0]:trans_d[locusfname][1]] #secret opt
                                print "s:", trans_d[locusfname][0], "e:", trans_d[locusfname][1], "trans end:", len(tempseq)
                            rhandle.seek(0)
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
                        seq.id = b[:-6].split("/")[-1][:-4]#[:5]
                        seq.name =""
                        seq.description =""
                        if opt == "-mss" and locusfname not in warninglist:
                            SeqIO.write(seq, fhandle, "fasta")
                            c1 += 1
                        elif opt == "-ms" or opt == "-me":
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
        for wle in warninglist:
            print wle
        totalloci += c1
    print "number of >90% close suboptimal hits", number
    print "number of loci with >90% close suboptimal hits", len(numberset)
    print "total number of loci found", totalloci
else:
    output = {} #main dctionary
    trans_d = {} #ahe loci and corresponding start and end of trans locus
    #reading blast file
    print "reading blastfile...", blastfilearg
    blastfile = open(blastfilearg, "rU")
    reader = csv.reader(blastfile, delimiter='\t')
    currentkey = ""
    best_eval = 0.0
    best_hit = 0
    # debug_file = open("debug.tab", "w")
    # debug_dict = {}
    # debug_c = 0
    for row in reader:
        if currentkey != row[0]: ##new query
            if float(row[10]) <= evalue:
                if row[0].split("/")[-1] in output: ##query is present - normaly never happens
                    print "warning: the key exists", row[0].split("/")[-1]
                else:
                    # debug_c = 0
                    output[row[0].split("/")[-1]] = row[1] ##standart output
                    #print "NORMAL"
                    best_eval = float(row[10])
                    best_hit = float(row[11])
                    # debug_c = 1
                    #output[AHE] = trans
                    #trans_d[AHE] = [start, finish] - for a given AHE, what was the start and end
                    if opt == "-s" or opt == "-ss":
                        start = min(int(row[8]), int(row[9])) #secret opt
                        end = max(int(row[8]), int(row[9])) #secret opt
                        trans_d[row[0].split("/")[-1]] = [start, end] #secret opt
                    print "--------------------------------------------------------"
                    print row[0].split("/")[-1], row[1], row[2], row[10], row[11]
        else: ##same query
            if currentmatch != row[1] and float(row[10]) <= evalue:
                if row[0].split("//")[-1] in output: 
                    #adopted for trinity assemblies
                    if row[1].split("_c")[0] == output[row[0].split("/")[-1]].split("_c")[0]:
                        print "warning: several isoforms detected", row[1], row[10], row[11]
                    else:
                        #print "warning: several matches detected", row[1], row[10], row[11], "; ratio with best_eval is", round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                        # print >> debug_file, row[0].split("/")[-1], row[1].split("_c")[0], math.log(float(row[10]))/math.log(best_eval)*100, float(row[11])/best_hit*100
                        # debug_c +=1
                        print "warning: several matches detected", row[1], row[10], row[11]#, "; ratio with best_eval is", round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                        if round(float(row[11])/best_hit*100) > 90: #round(math.log(float(row[10]))/math.log(best_eval)*100) > 90 or 
                            number += 1
                            numberset.add(row[0])
                else:
                    output[row[0].split("/")[-1]] = row[1] ##standart output
                    print "ADDITION"
        currentkey = row[0]
        currentmatch = row[1]
        currente = float(row[10])
        # debug_dict[currentkey.split("/")[-1]] = debug_c
    blastfile.close()
    # debug_file.close()
    count = int(len(output))
    print count, "targets found to be extracted"
    # debug_file2 = open("debug2.tab", "w")
    # for k, v in debug_dict.items():
    #     print >> debug_file2, k, v
    # debug_file2.close()
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
                    if opt == "-s" or opt == "-ss":
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
                    seq.id = trif.split("/")[-1][:-4]#[:5]
                    seq.name =""
                    seq.description =""
                    if opt == "-ss" and locusfname not in warninglist:
                        SeqIO.write(seq, fhandle, "fasta")
                        c1 += 1
                    elif opt == "-s":
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
    for wle in warninglist:
        print wle
print "done"