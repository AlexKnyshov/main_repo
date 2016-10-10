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
    print "FORMAT: python alexblparser.py [blastfile or folder] [asemblyfile or folder] [ahefolder] [evalue] [option: -n (normal), -s (extract only matched parts), -ss (short are discarded), -e (extended s option), -mn (multiple normal), -ms (multiple selected), -mss(short are discarded), -me (extended ms option)]"
    print "EXAMPLE: python alexblparser.py ./blast_outputs/ /transcriptomes/ ./fasta 1e-40 -mn"
    print "EXAMPLE: python alexblparser.py blast.tab trinity.fas ./fasta/ 1e-40 -n"
    sys.exit()

dash = "--------------------------------------------------------"

#function for creating an output folder. old stuff will be deleted
def mkdirfunc():
    if not os.path.exists ("./modified/"):
        os.makedirs("./modified") #creating folder if necessary
    else:
        shutil.rmtree("./modified/") #removing old files
        os.makedirs("./modified")

#function for copying alignment files before changing them
#all alignments are copied regardless of whether they will be modified or not
def copyfunc():
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

#function for parsing a blast output file
#returns a dictionary like this:
#output2[AHE] = [transf, transr, transb, ahef, aher, aheb, trans]
#AHE - query name, transf - target start pos, transr - target end, transb - forward or reverse target direction
#ahef - query start, aher - query end, aheb - query direction, trans - target name
#each query is allowed to have only 1 target
def readblastfilefunc(b, debugfile):
    print "processing", b
    print >> debugfile, "processing blastfile", b
    output2 = {}
    blastfile = open(b, "rU")
    reader = csv.reader(blastfile, delimiter='\t')
    currentkey = ""
    best_eval = 0.0
    best_hit = 0
    for row in reader:
        if currentkey != row[0]: ##new query
            if float(row[10]) <= evalue:
                if row[0].split("/")[-1] in output2: ##query is present
                    print "warning: the key exists", row[0].split("/")[-1]
                    print >> debugfile, "warning: the key exists", row[0].split("/")[-1]
                else:
                    print dash
                    print >> debugfile, dash
                    best_eval = float(row[10])
                    best_hit = float(row[11])
                    #check trans
                    if int(row[8]) < int(row[9]):
                        transb = True
                    else:
                        transb = False
                    transf = int(row[8])
                    transr = int(row[9])
                    #check query
                    if int(row[6]) < int(row[7]):
                        aheb = True
                    else:
                        aheb = False
                    ahef = int(row[6])
                    aher = int(row[7])
                    #store the data:
                    output2[row[0].split("/")[-1]] = [transf, transr, transb, ahef, aher, aheb, row[1]]

                    print "query", row[0].split("/")[-1], ", direction:", aheb, "; target", row[1], ", direction: ", transb
                    print >> debugfile, "query", row[0].split("/")[-1], ", direction:", aheb, "; target", row[1], ", direction: ", transb
                    print "identity:", row[2], ", eval:", row[10], ", bitscore", row[11]
                    print >> debugfile, "identity:", row[2], ", eval:", row[10], ", bitscore", row[11]
        else: ##same query
            if currentmatch != row[1] and float(row[10]) <= evalue:
                if row[0].split("/")[-1] in output2: 
                    #adopted for trinity assemblies
                    if row[1].split("_c")[0] == output2[row[0].split("/")[-1]][6].split("_c")[0]:
                        print "warning: several isoforms detected", row[1], row[10], row[11]
                        print >> debugfile, "warning: several isoforms detected", row[1], row[10], row[11]
                    else:
                        print "warning: several matches detected", row[1], row[10], row[11], "hit is", round(float(row[11])/best_hit*100), "%"#"; ratio with best_eval is"#, round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                        print >> debugfile, "warning: several matches detected", row[1], row[10], row[11], "hit is", round(float(row[11])/best_hit*100), "%"#"; ratio with best_eval is"#, round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                        if round(float(row[11])/best_hit*100) > 90:#round(math.log(float(row[10]))/math.log(best_eval)*100) > 90 or :
                            #number += 1
                            #numberset.add(row[0])
                            print >> debugfile, "number += 1"
                else:
                    #output[row[0].split("/")[-1]] = row[1] ##standart output
                    output2[row[0].split("/")[-1]] = [transf, transr, transb, ahef, aher, aheb, row[1]]
                    print >> debugfile, "ADDITION"
        currentkey = row[0]
        currentmatch = row[1]
        currente = float(row[10])
    blastfile.close()
    return output2

#function for determining the offsets of blast hit in a query and a target
#returns the list of four elements:
#gap - query hit start offset, gaprev - query hit end offset
#gapt - target hit start offset, gaptrev - target hit end offset
def gapfunc(output, locusfname, recs):
    if output[locusfname][5]:#aheb: #ahe is forward
        #print "AHE forw"
        gap = output[locusfname][3]#ahef #ahe start gap
        if output[locusfname][2]:#transb: #trans is forward
            ##need to check if trans has start - OK
            if output[locusfname][0] - gap >= 0:
                gapt = output[locusfname][0] - gap#transf - gap #this is the corrected trans start
            else:
                gapt = 0 # otherwise start from start
        else:
            ##need to check if trans is long enough - OK
            if output[locusfname][0] + gap <= len(tempseq):
                gapt = output[locusfname][0] + gap#transf + gap #this is the corrected trans start
            else:
                gapt = len(tempseq) # otherwise start from the end
        #check the ahe end:
        gaprev = len(recs[0].seq) - output[locusfname][4]#ali.get_alignment_length() - output[locusfname][0][4]
        #if gaprev >= 0:#aher #some AHE left
        if output[locusfname][2]:#transb: #trans is forward
            #check the trans end
            if output[locusfname][1] + gaprev <= len(tempseq): #check that trans is long enough
                gaptrev = output[locusfname][1] + gaprev#transf - gap #this is the corrected trans start
            else:
                gaptrev = len(tempseq)
        else:
            ##need to check if trans has start
            if output[locusfname][1] - gaprev >= 0:
                gaptrev = output[locusfname][1] - gaprev#transf + gap #this is the corrected trans start
            else:
                gaptrev = 0
    #reverse AHE:
    else: #ahe is reverse
        #print "AHE rev"
        gap = len(recs[0].seq) - output[locusfname][3]#ali.get_alignment_length() - output[locusfname][0][3]#ahef #ahe start gap
        if output[locusfname][2]:#transb: #trans is forward
            ##need to check if trans has start - OK
            if output[locusfname][0] - gap >= 0:
                gapt = output[locusfname][0] - gap#transf - gap #this is the corrected trans start
            else:
                gapt = 0
        else:
            ##need to check if trans is long enough - OK
            if output[locusfname][0] + gap <= len(tempseq):
                gapt = output[locusfname][0] + gap#transf + gap #this is the corrected trans start
            else:
                gapt = len(tempseq)
        #check the ahe end:
        gaprev = output[locusfname][4] # ahe end gap
        #if gaprev >= 0:#aher #some AHE left
        if output[locusfname][2]:#transb: #trans is forward
            #check the trans end
            if output[locusfname][1] + gaprev <= len(tempseq): #check that trans is long enough
                gaptrev = output[locusfname][1] + gaprev#transf - gap #this is the corrected trans start
            else:
                gaptrev = len(tempseq)
        else:
            ##need to check if trans has start
            if output[locusfname][1] - gaprev >= 0:
                gaptrev = output[locusfname][1] - gaprev#transf + gap #this is the corrected trans start
            else:
                gaptrev = 0
    return [gap, gaprev, gapt, gaptrev]

#function for preparing found target sequence for appending
#returns prepared sequence
def seqprepfunc(output, locusfname, opt, seq):
    rhandle = open("./modified/"+locusfname, "r")
    ali = SeqIO.parse(rhandle, "fasta")
    recs = list(ali)
    # identify gaps
    gaps = gapfunc(output, locusfname, recs)
    # print out trans and AHE stats...
    #print "BLAST:", "transf:",output[locusfname][0], "transr:",output[locusfname][1], "transb:",output[locusfname][2], "ahef:",output[locusfname][3], "aher:",output[locusfname][4], "aheb:",output[locusfname][5]
    #print "LOCUS", "ahe:", len(recs[0].seq), "trans:", len(seq.seq)
    #print gaps
    if opt == "-me" or opt == "-e":
        #TEST
        if output[locusfname][2]:
            #print "STATS:", output[locusfname][2], "start:", gaps[2],"stop:",  gaps[3]
            seq.seq = seq.seq[gaps[2]:gaps[3]]
        else:
            #print "STATS:", output[locusfname][2], "start:", gaps[3],"stop:",  gaps[2]
            seq.seq = seq.seq[gaps[3]:gaps[2]]
            seq = seq.reverse_complement()
    elif opt == "-ms" or opt == "-mss" or opt == "-s" or opt == "-ss":
        if output[locusfname][2] and output[locusfname][5]:
            seq.seq = seq.seq[output[locusfname][0]:output[locusfname][1]] #secret opt
        else:
            seq.seq = seq.seq[output[locusfname][1]:output[locusfname][0]]
            seq = seq.reverse_complement()
    elif opt == "-mn" or opt == "-n":
        if not output[locusfname][2] and not output[locusfname][5]:
            seq = seq.reverse_complement()
    rhandle.close()
    return seq

#function for appending the seqeunce to the alignemnt
#returns 1 if success, 0 otherwise
def seqwritefunc(seq, opt, locusfname, seqname):
    fhandle = open("./modified/"+locusfname, "a")
    seq.id = seqname
    seq.name =""
    seq.description =""
    if (opt == "-mss" or opt == "-ss") and (float(len(seq.seq)) / loclenfunc(locusfname) > 0.8):
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt == "-ms" or opt == "-me" or opt == "-s" or opt == "-e":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt == "-mn" or opt == "-n":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    else:
        return 0
    fhandle.close()

#function for determining the length of the query alignment
#returns length (int)
def loclenfunc(locusfname):
    rhandle = open("./modified/"+locusfname, "r")
    ali = list(SeqIO.parse(rhandle, "fasta"))
    return len(ali[0])
    rhandle.close()
#---------------------------------------------------------------------

print "alexblparser run with option", opt, "selected"
debugfile = open("alexblparser.log", "w")
print >> debugfile, "debug file start"
print >> debugfile, "command line parameters:", sys.argv
#make modified dir
print >> debugfile, "make modified dir..."
mkdirfunc()

#copy files
print >> debugfile, "copy files..."
copyfunc()

#multi db option
if opt == "-mn" or opt == "-ms" or opt == "-mss" or opt == "-me":
    #reading the blastfile
    blastlist = glob.glob(blastfilearg+"/*.blast")
    translist = glob.glob(trif+"/*.fasta")
else:
    blastlist = [blastfilearg]
    translist = [trif]

print "list of transcriptomes:"
print >> debugfile, "list of transcriptomes:"
for l in translist:
    print l
    print >> debugfile, l
#debug vars
number = 0
numberset = set()
totalloci = 0
#parsing blast files
print "parsing blast files..."
print >> debugfile, "parsing blast files..."
for b in blastlist:
    #obtain a dictionary with blast results
    output = readblastfilefunc(b, debugfile)
    print dash
    print >> debugfile, dash
    count = int(len(output))
    print count, "targets found to be extracted"
    print >> debugfile, count, "targets found to be extracted"
    
    print "scanning the transcriptome..."
    print >> debugfile, "scanning the transcriptome..."
    warninglist = []
    #get the transcriptome filename, matching blast filename
    for t_file in translist:
        if b[:-6].split("/")[-1] in t_file:
            #print >> debugfile, b[:-6].split("/")[-1], t_file
            inputf = SeqIO.parse(t_file, "fasta")
            seqname = b[:-6].split("/")[-1]
            transname = t_file.split("/")[-1]
            break
    print >> debugfile, "target:", transname, "; target name:", seqname
    if not inputf:
        print "error, transcriptome file is not found"
        print >> debugfile, "error, transcriptome file is not found"
        break
    print "searching for contigs in:", transname
    print >> debugfile, "searching for contigs in:", transname
    c1 = 0
    extr_loci = []                                                  
    for seq in inputf:
        #print count
        if count == 0:
            print dash
            print >> debugfile, dash
            print "search terminated"
            print >> debugfile, "search terminated"
            break
        else:
            # for val in output.values():
            #     if seq.id in val: #if contig is in ahe (was found as a blast hit)
            #         print dash
            #         print >> debugfile, dash
            #         print "target", seq.id, "found as a blast hit, looking up the locus..."
            #         print >> debugfile, "target", seq.id, "found as a blast hit, looking up the locus..."
            for locusfname,y in output.items(): #checking all ahe (seach which ahe has it)
                #locusfname - AHE name
                #y - trans locus name
                #seq.id = trans locus name as well
                if y[6] == seq.id and locusfname not in extr_loci:
                    print dash
                    print >> debugfile, dash
                    print "target", seq.id, "found as a blast hit, matched the locus", locusfname
                    print >> debugfile, "target", seq.id, "found as a blast hit, matched the locus",locusfname
                    #print "locus", locusfname
                    #print >> debugfile, "locus", locusfname
                    temp = seq.id
                    tempseq = seq.seq #secret opt
                    
                    #check direction and length
                    seq = seqprepfunc(output, locusfname, opt, seq)
                    print "extracted length:", len(seq.seq)
                    print >> debugfile, "extracted length:", len(seq.seq)
                    #append sequence
                    if seqwritefunc(seq, opt, locusfname, seqname) == 1:
                        c1 += 1
                    elif opt == "-mss" or opt == "-ss":
                        print "Warning: blast hit is too short"
                        print >> debugfile, "Warning: blast hit is too short"
                        warninglist.append(locusfname)
                    extr_loci.append(locusfname)
                    seq.id = temp
                    seq.seq = tempseq #secret opt
                    count -= 1
                # if y[6] == seq.id and locusfname in extr_loci:
                #     print >> debugfile, "locus is already extracted"
                #     print "locus is already extracted"
    count = int(len(output))
    print "queries found:", len(extr_loci)
    print >> debugfile, "queries found:", len(extr_loci)
    print c1, "sequences extracted"
    print >> debugfile, c1, "seqeunces extracted"
    print "warning list:", len(warninglist)
    print >> debugfile, "warning list:", len(warninglist)
    for wle in warninglist:
        print wle
        print >> debugfile, wle
    totalloci += c1

print "number of >90% close suboptimal hits", number
print "number of queries with >90% close suboptimal hits", len(numberset)
print "total number of queries extracted", totalloci, ", in average", totalloci/float(len(translist)), "per database"

print >> debugfile, "number of >90% close suboptimal hits", number
print >> debugfile, "number of queries with >90% close suboptimal hits", len(numberset)
print >> debugfile, "total number of queries extracted", totalloci, ", in average", totalloci/float(len(translist)), "per database"
print >> debugfile, "done"
debugfile.close()
print "done"