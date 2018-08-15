from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os
import shutil
import csv
import sys
import math
if len(sys.argv) == 7:
    blastfilearg = sys.argv[1]
    trif = sys.argv[2]
    ahefoldarg = sys.argv[3]
    evalue = float(sys.argv[4])
    opt = sys.argv[5]
    paropt = int(sys.argv[6])
    if paropt < 1:
        print "[paropt] cannot be less than 1"
        sys.exit()
else:
    print "FORMAT: python mainblparser.py [blastfile or folder] [asemblyfile or folder] [ahefolder] [evalue] [option: -n (normal), -s (extract only matched parts), -ss (short are discarded), -e[value] (extended s option), -a (blast region extraction and stiching), -b (extract between two outmost blast regions), -mn (multiple normal), -ms (multiple selected), -mss(short are discarded), -me[value] (extended ms option), -ma, -mb] [number of hits extracted per query]"
    print "EXAMPLE: python mainblparser.py ./blast_outputs/ /transcriptomes/ ./fasta 1e-40 -mn 2"
    print "EXAMPLE: python mainblparser.py blast.tab trinity.fas ./fasta/ 1e-40 -n 1"
    print "HELP for -e / -me option: [value] indicates a size of the flank in bp that needs to be extracted"
    print "HELP for -e / -me option: when [value]=0, extracts flanks to match the length of the query"
    print "EXAMPLE: python mainblparser.py blast.tab trinity.fas ./fasta/ 1e-40 -e1500 1"
    print "1500 bp from each side (if possible) will be added to blast hit region and extracted"
    sys.exit()

dash = "--------------------------------------------------------"
starterr = "Start offset is too large"
enderr = "End offset is too large"
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
    for x in glob.glob(ahefoldarg+"/*.fas*"):
        locusfname = x.split("/")[-1]
        #print locusfname
        if not os.path.exists ("./modified/"+locusfname):
            prog = "copying "+str(locusfname)+"..."
            sys.stdout.write(prog+"\r")
            sys.stdout.flush()
            shutil.copy2(ahefoldarg+locusfname, "./modified")
    print ""

#function for parsing a blast output file
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
    query_counter = 0
    for row in reader:
        if float(row[10]) <= evalue:
            if currentkey != row[0]: ##new query
                print >> debugfile, "### NEW QUERY & TARGET###"
                query_counter = 0
                parts_counter = 0
                #if float(row[10]) <= evalue:
                if row[0].split("/")[-1] in output2: ##query is present
                    print "warning: the key exists", row[0].split("/")[-1]
                    print >> debugfile, "warning: the key exists", row[0].split("/")[-1]
                else:
                    print dash
                    print >> debugfile, dash
                    best_eval = float(row[10])
                    best_hit = float(row[11])
                    output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)] = rowfunc(row)
                    print "query", row[0].split("/")[-1], ", direction:", output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)][5], "; target", row[1], ", direction: ", output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)][2]
                    print >> debugfile, "query", row[0].split("/")[-1], ", direction:", output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)][5], "; target", row[1], ", direction: ", output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)][2]
                    print "identity:", row[2], ", eval:", row[10], ", bitscore", row[11]
                    print >> debugfile, "identity:", row[2], ", eval:", row[10], ", bitscore", row[11]
                    #query_counter += 1
            else: ##same query
                print >> debugfile, "### SAME QUERY ###"
                if currentmatch != row[1]:# and float(row[10]) <= evalue:#new target
                    print >> debugfile, "### NEW TARGET ###"
                    query_counter += 1
                    parts_counter = 0
                    if row[0].split("/")[-1]+"-copy0"+"-copy"+str(parts_counter) in output2: #
                        #adopted for trinity assembly sequence headers
                        if row[1].split("_c")[0] == output2[row[0].split("/")[-1]+"-copy0"+"-copy"+str(parts_counter)][6].split("_c")[0]:
                            print "warning: several isoforms detected", row[1], row[10], row[11]
                            print >> debugfile, "warning: several isoforms detected", row[1], row[10], row[11]
                        else:
                            print "warning: several matches detected", row[1], row[10], row[11], "hit is", round(float(row[11])/best_hit*100), "%"#"; ratio with best_eval is"#, round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                            print >> debugfile, "warning: several matches detected", row[1], row[10], row[11], "hit is", round(float(row[11])/best_hit*100), "%"#"; ratio with best_eval is"#, round(math.log(float(row[10]))/math.log(best_eval)*100), "% (log), hit is", round(float(row[11])/best_hit*100), "%"
                            if query_counter < paropt:#round(float(row[11])/best_hit*100) > 70 and query_counter < paropt:#round(math.log(float(row[10]))/math.log(best_eval)*100) > 90 or :
                                print >> debugfile, "homolog added"
                                #store the additional data:
                                #query_counter += 1
                                output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)] = rowfunc(row)
                                print "locus", row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)
                            else:
                                print >> debugfile, "Too many copies, ignored", query_counter, parts_counter, paropt
                        #query_counter += 1
                    else:
                        output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)] = rowfunc(row)
                        print >> debugfile, "ERROR: if this appears, the script may have worked incorrectly and/or the input file is formatted incorrectly"
                else:
                    #if float(row[10]) <= evalue:
                    if query_counter < paropt:                    
                        print >> debugfile, "same query and target, current interval", row[1], int(row[8]), int(row[9]),"key", row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter), "stored interval", output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy0"][0:2]
                        new = True
                        for key,value in output2.items():#create comparative dict
                            comp_dict = {}
                            num = 0
                            if "-copy".join(key.split("-copy")[0:2]) == row[0].split("/")[-1]+"-copy"+str(query_counter) and value[6] == row[1]:
                                if getOverlap([int(row[8]), int(row[9])],value[0:2]) > 0 or getOverlap([int(row[6]), int(row[7])],value[3:5]) > 0:
                                    new = False
                                    #correct region to elongate hit region:

                        if new:
                            parts_counter += 1
                            print >> debugfile, "interval equal 0, new region recovered, value7", parts_counter
                            output2[row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter)] = rowfunc(row)
                            print >> debugfile, row[0].split("/")[-1]+"-copy"+str(query_counter)+"-copy"+str(parts_counter), rowfunc(row)
                        else:
                            print >> debugfile, "interval greater than 0, rowfunc not called, value7", parts_counter
                    else:
                        print >> debugfile, "(additional stiching) Too many copies, ignored", query_counter, parts_counter, paropt
            currentkey = row[0]
            currentmatch = row[1]
            currente = float(row[10])
        else:
            print >> debugfile, "low evalue, skipping row", row
    blastfile.close()
    return output2

#function to return a list for dict
#returns a dictionary like this:
#output2[AHE] = [transf, transr, transb, ahef, aher, aheb, row[1]]
#AHE - query name, transf - target start pos, transr - target end, transb - forward or reverse target direction
#ahef - query start, aher - query end, aheb - query direction, trans - target name
def rowfunc(row):
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
    return [transf, transr, transb, ahef, aher, aheb, row[1]]

#function to compute overlap between two ranges supplied as lists with start and end
#returns overlap value
def getOverlap(a, b):
    a0=min(a)
    a1=max(a)
    b0=min(b)
    b1=max(b)
    #print >> debugfile, "interval computation", a, b, a0, a1, b0, b1, max(0, min(a1, b1) - max(a0, b0))
    return max(0, min(a1, b1) - max(a0, b0))

#function to count amount of blast hits for each copy
#returns dict[locus] = copy_number
def partsnumfunc(input_dict):
    output_dict = {}
    for key1, value1 in input_dict.items():
        if "-copy".join(key1.split("-copy")[0:2]) in output_dict:
            output_dict["-copy".join(key1.split("-copy")[0:2])] += 1
        else:
            output_dict["-copy".join(key1.split("-copy")[0:2])] = 1
    return output_dict


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
                print starterr
                gapt = 0 # otherwise start from start
        else:
            ##need to check if trans is long enough - OK
            if output[locusfname][0] + gap <= len(tempseq):
                gapt = output[locusfname][0] + gap#transf + gap #this is the corrected trans start
            else:
                print enderr
                gapt = len(tempseq) # otherwise start from the end
        #check the ahe end:
        gaprev = len(recs[0].seq) - output[locusfname][4]#ali.get_alignment_length() - output[locusfname][0][4]
        #if gaprev >= 0:#aher #some AHE left
        if output[locusfname][2]:#transb: #trans is forward
            #check the trans end
            if output[locusfname][1] + gaprev <= len(tempseq): #check that trans is long enough
                gaptrev = output[locusfname][1] + gaprev#transf - gap #this is the corrected trans start
            else:
                print enderr
                gaptrev = len(tempseq)
        else:
            ##need to check if trans has start
            if output[locusfname][1] - gaprev >= 0:
                gaptrev = output[locusfname][1] - gaprev#transf + gap #this is the corrected trans start
            else:
                print starterr
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
                print starterr
                gapt = 0
        else:
            ##need to check if trans is long enough - OK
            if output[locusfname][0] + gap <= len(tempseq):
                gapt = output[locusfname][0] + gap#transf + gap #this is the corrected trans start
            else:
                print enderr
                gapt = len(tempseq)
        #check the ahe end:
        gaprev = output[locusfname][4] # ahe end gap
        #if gaprev >= 0:#aher #some AHE left
        if output[locusfname][2]:#transb: #trans is forward
            #check the trans end
            if output[locusfname][1] + gaprev <= len(tempseq): #check that trans is long enough
                gaptrev = output[locusfname][1] + gaprev#transf - gap #this is the corrected trans start
            else:
                print enderr
                gaptrev = len(tempseq)
        else:
            ##need to check if trans has start
            if output[locusfname][1] - gaprev >= 0:
                gaptrev = output[locusfname][1] - gaprev#transf + gap #this is the corrected trans start
            else:
                print starterr
                gaptrev = 0
    return [gap, gaprev, gapt, gaptrev]



#gap function for -e and -me options
def gapfunc_e(output, locusfname, e_range):
    if output[locusfname][5]:#aheb: #ahe is forward
        #print "AHE forw"
        if output[locusfname][2]:#transb: #trans is forward
            ##need to check if trans has start - OK
            if output[locusfname][0] - e_range >= 0:
                gapt = output[locusfname][0] - e_range#transf - gap #this is the corrected trans start
            else:
                print starterr
                gapt = 0 # otherwise start from start
        else:
            ##need to check if trans is long enough - OK
            if output[locusfname][0] + e_range <= len(tempseq):
                gapt = output[locusfname][0] + e_range#transf + gap #this is the corrected trans start
            else:
                gapt = len(tempseq) # otherwise start from the end
                print enderr
        #check the ahe end:
        #if gaprev >= 0:#aher #some AHE left
        if output[locusfname][2]:#transb: #trans is forward
            #check the trans end
            if output[locusfname][1] + e_range <= len(tempseq): #check that trans is long enough
                gaptrev = output[locusfname][1] + e_range#transf - gap #this is the corrected trans start
            else:
                print enderr
                gaptrev = len(tempseq)
        else:
            ##need to check if trans has start
            if output[locusfname][1] - e_range >= 0:
                gaptrev = output[locusfname][1] - e_range#transf + gap #this is the corrected trans start
            else:
                print starterr
                gaptrev = 0
    #reverse AHE:
    else: #ahe is reverse
        #print "AHE rev"
        if output[locusfname][2]:#transb: #trans is forward
            ##need to check if trans has start - OK
            if output[locusfname][0] - e_range >= 0:
                gapt = output[locusfname][0] - e_range#transf - gap #this is the corrected trans start
            else:
                gapt = 0
                print starterr
        else:
            ##need to check if trans is long enough - OK
            if output[locusfname][0] + e_range <= len(tempseq):
                gapt = output[locusfname][0] + e_range#transf + gap #this is the corrected trans start
            else:
                gapt = len(tempseq)
                print enderr
        #check the ahe end:
        #if gaprev >= 0:#aher #some AHE left
        if output[locusfname][2]:#transb: #trans is forward
            #check the trans end
            if output[locusfname][1] + e_range <= len(tempseq): #check that trans is long enough
                gaptrev = output[locusfname][1] + e_range#transf - gap #this is the corrected trans start
            else:
                gaptrev = len(tempseq)
                print enderr
        else:
            ##need to check if trans has start
            if output[locusfname][1] - e_range >= 0:
                gaptrev = output[locusfname][1] - e_range#transf + gap #this is the corrected trans start
            else:
                print starterr
                gaptrev = 0
    return [gapt, gaptrev]




#function for preparing found target sequence for appending
#returns prepared sequence
def seqprepfunc(output, locusfname, opt, seq):
    if opt[:3] == "-me" or opt[:2] == "-e":
        if opt[:3] == "-me" and int(opt[3:]) == 0  or opt[:2] == "-e" and int(opt[2:]) == 0 :
            rhandle = open("./modified/"+locusfname.split("-copy")[0], "r")
            ali = SeqIO.parse(rhandle, "fasta")
            recs = list(ali)
            # identify gaps
            gaps = gapfunc(output, locusfname, recs)
            if output[locusfname][2] and output[locusfname][5]:
                seq.seq = seq.seq[gaps[2]:gaps[3]]
            elif not output[locusfname][2] and not output[locusfname][5]:
                seq.seq = seq.seq[gaps[3]:gaps[2]]
            else:
                if gaps[2] < gaps[3]:
                    seq.seq = seq.seq[gaps[2]:gaps[3]]
                else:
                    seq.seq = seq.seq[gaps[3]:gaps[2]]
                seq = seq.reverse_complement()
            rhandle.close()
        else:
            # identify gaps
            if opt[:3] == "-me":
                gapsE = gapfunc_e(output, locusfname, int(opt[3:]))
            else: #-e option
                gapsE = gapfunc_e(output, locusfname, int(opt[2:]))
            if output[locusfname][2] and output[locusfname][5]:
                seq.seq = seq.seq[gapsE[0]:gapsE[1]]
            elif not output[locusfname][2] and not output[locusfname][5]:
                seq.seq = seq.seq[gapsE[1]:gapsE[0]]
            else:
                if gapsE[0] < gapsE[1]:
                    seq.seq = seq.seq[gapsE[0]:gapsE[1]]
                else:
                    seq.seq = seq.seq[gapsE[1]:gapsE[0]]
                seq = seq.reverse_complement()
    elif opt == "-ms" or opt == "-mss" or opt == "-s" or opt == "-ss" or opt == "-msl" or opt == "-sl":
        if output[locusfname][2] and output[locusfname][5]:
            seq.seq = seq.seq[output[locusfname][0]:output[locusfname][1]]
        elif not output[locusfname][2] and not output[locusfname][5]:
            seq.seq = seq.seq[output[locusfname][1]:output[locusfname][0]]
        else:
            if output[locusfname][0] > output[locusfname][1]:
                seq.seq = seq.seq[output[locusfname][1]:output[locusfname][0]]
            else:
                seq.seq = seq.seq[output[locusfname][0]:output[locusfname][1]]
            seq = seq.reverse_complement()
    elif opt == "-mn" or opt == "-n":
        if (output[locusfname][2] and not output[locusfname][5]) or (not output[locusfname][2] and output[locusfname][5]):
            seq = seq.reverse_complement()
    return seq

#alternative function for preparing target sequence for appending
#stiches non-overlapping blast hit regions
#returns prepared sequence
def seqprepfunc_es(output, locus1, opt, seq1, partsnumdict):
    if opt == "-ma" or opt == "-a":
        rhandle = open("./modified/"+locusfname.split("-copy")[0], "r")
        ali = SeqIO.parse(rhandle, "fasta")
        recs = list(ali)
        startlist = {}
        endlist = {}
        outseq = Seq("")
        endpos = 0
        indel = 0
        print >> debugfile, "locus:", locusfname, "num parts:", partsnumdict[locus1]
        for loc in output.keys():
            if "-copy".join(loc.split("-copy")[0:2]) == locus1:
                if output[loc][5]:
                    startlist[loc] = output[loc][3]-1
                    endlist[loc] = output[loc][4]
                else:
                    startlist[loc] = output[loc][4]-1
                    endlist[loc] = output[loc][3]
        #print >> debugfile, "DEBUG1", locus1, "startlist and endlist", startlist, endlist
        for key1 in sorted(startlist, key=startlist.__getitem__):

            if endpos == startlist[key1]:
                print >> debugfile, "DEBUG1", key1, "no start added gap", endpos, startlist[key1]
            else:
                outseq += Seq("-"*(startlist[key1]-endpos))
                #print >> debugfile, "DEBUG1", key1, "start added gap =", endpos, startlist[key1], startlist[key1]-endpos

            print >> debugfile, "locus:", startlist[key1], "output", output[key1]
            if output[key1][2] and output[key1][5]:
                outseq += seq1[(output[key1][0]-1):(output[key1][1])]
            elif not output[key1][2] and not output[key1][5]:
                outseq += seq1[(output[key1][1]-1):(output[key1][0])]
            else:
                if output[key1][0] > output[key1][1]:#transrev
                    outseq += seq1[(output[key1][1]-1):(output[key1][0])].reverse_complement()
                else:#transforw
                    outseq += seq1[(output[key1][0]-1):(output[key1][1])].reverse_complement()
                #outseq = outseq.reverse_complement()
            endpos = endlist[key1]
            #outseq += Seq("NNNNN")
            #print >> debugfile, "DEBUG1", key1, "startpos and endpos", endpos
        outseq += Seq("-"*(len(recs[0].seq)-endpos))
    elif opt == "-mb" or opt == "-b":
        startlist = {}
        endlist = {}
        #outseq = Seq("")
        print >> debugfile, "locus:", locusfname, "num parts:", partsnumdict[locus1]
        for loc in output.keys():
            if "-copy".join(loc.split("-copy")[0:2]) == locus1:
                if output[loc][2]:
                    startlist[loc] = output[loc][0]
                    endlist[loc] = output[loc][1]
                else:
                    startlist[loc] = output[loc][1]
                    endlist[loc] = output[loc][0]
        startpos = min(startlist.values())
        endpos = max(endlist.values())
        print >> debugfile, "locus:", locus1, "coords", startpos, endpos
        if output[locus1+"-copy0"][2] and output[locus1+"-copy0"][5]:
            outseq = seq1[startpos:endpos]
        elif not output[locus1+"-copy0"][2] and not output[locus1+"-copy0"][5]:
            #outseq = seq1[endpos:startpos]
            outseq = seq1[startpos:endpos]
        else:
            outseq = seq1[startpos:endpos].reverse_complement()
            # if output[locus1+"-copy0"][0] > output[locus1+"-copy0"][1]:
            #     outseq = seq1[endpos:startpos].reverse_complement()
            # else:
            #     outseq = seq1[startpos:endpos].reverse_complement()
            #outseq = outseq.reverse_complement()
    return outseq#+Seq("NNNNN")

#function for appending the seqeunce to the alignemnt
#returns 1 if success, 0 otherwise
def seqwritefunc(seq, opt, locusfname, seqname, paropt):
    fhandle = open("./modified/"+locusfname.split("-copy")[0], "a")
    if paropt == 1:
        seq.id = seqname
    else:
        seq.id = seqname+"."+locusfname.split("-copy")[1]
    seq.name =""
    seq.description =""
    if (opt == "-mss" or opt == "-ss") and (float(len(seq.seq)) / loclenfunc(locusfname) > 0.8):
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt == "-ms" or opt == "-s":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt[:3] == "-me" or opt[:2] == "-e":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt == "-mn" or opt == "-n":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt == "-ma" or opt == "-a":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    elif opt == "-mb" or opt == "-b":
        SeqIO.write(seq, fhandle, "fasta")
        return 1
    else:
        return 0
    fhandle.close()

#alternative function for writing, all in 1 file
def altwritefunc(seq, opt, locusfname, seqname, paropt):
    fhandle = open("./"+seqname+"_conSeqs.fasta", "a")
    if paropt == 1:
        seq.id = locusfname.split("_")[-1].split(".")[0]
    else:
        seq.id = locusfname.split("_")[-1].split(".")[0]+"."+locusfname.split("-copy")[1]
    seq.name =""
    seq.description =""
    print >> fhandle, ">"+seq.id
    print >> fhandle, seq.seq
    return 1
    fhandle.close()

#function for determining the length of the query alignment
#returns length (int)
def loclenfunc(locusfname):
    rhandle = open("./modified/"+locusfname.split("-copy")[0], "r")
    ali = list(SeqIO.parse(rhandle, "fasta"))
    return len(ali[0])
    rhandle.close()
#---------------------------------------------------------------------

print "mainblparser run with option", opt, "selected"
debugfile = open("mainblparser.log", "w")
print >> debugfile, "debug file start"
print >> debugfile, "command line parameters:", sys.argv
#make modified dir
print >> debugfile, "make modified dir..."
mkdirfunc()

#copy files
print >> debugfile, "copy files..."
copyfunc()

#multi db option
if opt == "-mn" or opt == "-ms" or opt == "-mss" or opt[:3] == "-me" or opt == "-msl" or opt == "-ma" or opt == "-mb":
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
b1 = 0
for b in blastlist:
    b1 += 1
    #obtain a dictionary with blast results
    output = readblastfilefunc(b, debugfile)
    partsnumdict = partsnumfunc(output)
    print >> debugfile, partsnumdict
    print dash
    print >> debugfile, dash
    count = int(len(partsnumdict))
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
    if opt == "-msl" and os.path.exists("./"+seqname+"_conSeqs.fasta"):
        os.remove("./"+seqname+"_conSeqs.fasta")
    for seq in inputf:
        #print count
        if count == 0:
            print dash
            print >> debugfile, dash
            print "search terminated"
            print >> debugfile, "search terminated"
            break
        else:
            #for locusfname,y in output.items(): #checking all ahe (seach which ahe has it)
            for locus1 in partsnumdict.keys():
                #locusfname - AHE name
                #y - trans locus name
                #seq.id = trans locus name as well
                locusfname = locus1+"-copy0"
                y = output[locusfname]
                if y[6] == seq.id and locusfname not in extr_loci:
                    print dash
                    print >> debugfile, dash
                    print "target", seq.id, "found as a blast hit, matched the locus", locus1
                    print >> debugfile, "target", seq.id, "found as a blast hit, matched the locus",locus1
                    temp = seq.id
                    tempseq = seq.seq
                    temprec = seq
                    ######
                    print >> debugfile, "partsnumdict", partsnumdict[locus1]
                    if opt == "-a" or opt == "-b" or opt == "-ma" or opt == "-mb":
                        print "number of fragments to stich", partsnumdict[locus1]
                        print >> debugfile, "number of fragments to stich", partsnumdict[locus1]
                        writeseq = seqprepfunc_es(output, locus1, opt, temprec, partsnumdict)
                    else:
                        writeseq = seqprepfunc(output, locus1+"-copy0", opt, temprec)
                    extr_loci.append(locus1+"-copy0")
                    # if partsnumdict[locus1] == 1:
                    #     #print >> debugfile, "SEQ", seq.seq, len(seq.seq)
                    #     #print >> debugfile, "output debug",output[locus1+"-copy"+str(part)]#, seqprepfunc(output, locus1+"-copy"+str(part), opt, seq, partsnumdict)
                    #     if opt == "-a" or opt == "-b" or opt == "-ma" or opt == "-mb":
                    #     #writeseq = seqprepfunc(output, locus1+"-copy0", opt, temprec)
                    #         writeseq = seqprepfunc_es(output, locus1, opt, temprec, partsnumdict)
                    #     else:
                    #         writeseq = seqprepfunc(output, locus1+"-copy0", opt, temprec)
                    #     #else:
                    #     #    writeseq = writeseq + Seq("N"*5) + seqprepfunc(output, locus1+"-copy"+str(part), opt, temprec, partsnumdict)
                    #     #print >> debugfile, "output min", min(output[locus1+"-copy"+str(part)][0:2])
                    #     extr_loci.append(locus1+"-copy0")
                    #     #count -= 1
                    #     #print >> debugfile, writeseq, len(writeseq)
                    # else:
                    #     #print >> debugfile, "SEQ", seq.seq, len(seq.seq)
                    #     #print >> debugfile, "output debug",output[locus1+"-copy"+str(part)]#, seqprepfunc(output, locus1+"-copy"+str(part), opt, seq, partsnumdict)
                    #     #if part == 0:
                    #     writeseq = seqprepfunc_es(output, locus1, opt, temprec, partsnumdict)
                    #     #else:
                    #     #    writeseq = writeseq + Seq("N"*5) + seqprepfunc(output, locus1+"-copy"+str(part), opt, temprec, partsnumdict)
                    #     #print >> debugfile, "output min", min(output[locus1+"-copy"+str(part)][0:2])
                    #     extr_loci.append(locus1+"-copy0")
                    #     #count -= partsnumdict[locus1]
                    #     #print >> debugfile, writeseq, len(writeseq)

                    #check direction and length
                    #seq = seqprepfunc(output, locusfname, opt, seq, partsnumdict)
                    print "extracted length:", len(str(writeseq.seq).replace("-", ""))
                    print >> debugfile, "extracted length:", len(str(writeseq.seq).replace("-", ""))
                    #append sequence
                    if opt == "-msl" or opt == "-sl":
                        altwritefunc(writeseq, opt, locusfname, seqname, paropt)
                        c1 += 1
                    else:
                        if seqwritefunc(writeseq, opt, locusfname, seqname, paropt) == 1:
                           c1 += 1
                        elif opt == "-mss" or opt == "-ss":
                            print "Warning: blast hit is too short"
                            print >> debugfile, "Warning: blast hit is too short"
                            warninglist.append(locusfname)
                    #extr_loci.append(locusfname)
                    seq.id = temp
                    seq.seq = tempseq
                    count -= 1
                    print "progress:", count, "loci left, working on", transname, b1, "/", len(blastlist)
    if opt == "-msl" or opt == "-sl":
        print "sorting file....."
        fhandle = open("./"+seqname+"_conSeqs.fasta", "r")
        transdict = {}
        for seq in SeqIO.parse(fhandle, "fasta"):
            transdict[float(seq.id[1:])] = seq.seq
        fhandle.close()
        fhandle2 = open("./"+seqname+"_conSeqs.fasta", "w")
        for key, value in sorted(transdict.items()):
            print >> fhandle2, ">L"+str(key)
            print >> fhandle2, value
        fhandle2.close()
    count = int(len(partsnumdict))
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