from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import os
import shutil
import csv
import sys
import math
import itertools
import datetime
if len(sys.argv) == 10:
    blastfilearg = sys.argv[1]
    targetf = sys.argv[2]
    queryf = sys.argv[3] #add possibility to interpret this as query file, split it per locus and 
    evalue = float(sys.argv[4])
    filefolder = sys.argv[5]
    extractiontype = sys.argv[6]
    if extractiontype == "-e":
        print "option -e must have a numerical value, for ex. -e100 would extract 100bp flanks"
        sys.exit()
    contignum = int(sys.argv[7])
    if contignum == 0:
        print "[contignum] set to 0, extracting all contigs"
    if sys.argv[8] == "-Ry":
        reciprocate = True
    else:
        reciprocate = False
    if sys.argv[9] == "-IMy":
        interstich = True
    else:
        interstich = False
    print blastfilearg, targetf, queryf, evalue, "filefolder:",filefolder, "extractiontype:",extractiontype, "contignum:",contignum, "reciprocate:", reciprocate, "interstich:",interstich
else:
    print "FORMAT: python mainblparser.py [blast file or folder] [target file or folder] [query folder] [evalue] [single (-S)/ multiple file mode (-M)] [extraction type] [number of target contigs per query (if 0, extract all)] [check reciprocal best match] [perform intercontig stiching]"
    print ""
    print "Extraction types: -n (normal), -s (only best hit region), -e[value] (only best hit region plus flanks in bp), -a (extract all hit regions and join them), -b (extract region between two outmost blast regions)"
    print ""
    print "EXAMPLE: python mainblparser.py ./blast_outputs/ /transcriptomes/ ./fasta 1e-40 -M -n 2 -Ry -IMy"
    print "EXAMPLE: python mainblparser.py blast.tab trinity.fas ./fasta/ 1e-40 -S -n 1 -Rn -IMn"
    sys.exit()

#reciprocate = True
dash = "--------------------------------------------------------"
# starterr = "Start offset is too large"
# enderr = "End offset is too large"
warninglist = []

def messagefunc(msg, f, fl=True):
    if fl:
        sys.stdout.write(msg+"\r")
        sys.stdout.flush()
    else:
        print ""
        print msg
    print >> f, msg


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
    messagefunc("copying files *.fa*", debugfile, False)
    copyfunc_c = 0
    for x in glob.glob(queryf+"/*.fa*"):
        locusfname = x.split("/")[-1]
        #print locusfname
        if not os.path.exists ("./modified/"+locusfname):
            prog = "copying "+str(locusfname)+"..."
            messagefunc(prog, debugfile)
            shutil.copy2(queryf+"/"+locusfname, "./modified")
            copyfunc_c += 1
    messagefunc("copied "+str(copyfunc_c)+" files", debugfile, False)

#function for parsing a blast output file
#each query is allowed to have only 1 target
def readblastfilefunc(b, debugfile):
    messagefunc("processing "+b, debugfile, False)
    querydict = {}
    targetdict = {}
    blastfile = open(b, "rU")
    reader = csv.reader(blastfile, delimiter='\t')
    linecounter = 0
    for row in reader:
        if float(row[10]) <= evalue:
            qname = row[0].split("/")[-1]
            tname = row[1]
            #populate query table
            if qname in querydict:
                if tname in querydict[qname]:
                    querydict[qname][tname][linecounter] = rowfunc(row)
                else:
                    querydict[qname][tname] = {linecounter: rowfunc(row)}
            else:
                querydict[qname] = {tname: {linecounter: rowfunc(row)}}
            #populate target table
            if tname in targetdict:
                if qname in targetdict[tname]:
                    targetdict[tname][qname][linecounter] = rowfunc(row)
                else:
                    targetdict[tname][qname] = {linecounter: rowfunc(row)}
            else:
                targetdict[tname] = {qname: {linecounter: rowfunc(row)}}

        # else:
        #     print >> debugfile, "low evalue, skipping row", row
        linecounter += 1
    blastfile.close()
    return querydict, targetdict


#implement reciprocator break value, default 10. DO percet, e.g. 0.1 of the shortes range
def reciprocator(inpdict, query, range1, range2, emax, bitscore, target):
    messagefunc("running reciprocator on target "+target, debugfile)
    cond = True
    for key, val in inpdict.items():
        if key != query: #all other queries
            target_ranks_temp = compute_ranks(val) #get all pieces, compute values for this query
            if getOverlap([target_ranks_temp[2],target_ranks_temp[3]],[range1,range2]) > 10:
                if target_ranks_temp[0] < emax:
                    #print key, target_ranks_temp[2],target_ranks_temp[3], target_ranks_temp[0], query, range1,range2, emax
                    cond = False
                    break
                elif target_ranks_temp[0] == emax:
                    if target_ranks_temp[1] > bitscore:
                        cond = False
                        break
                    else:
                        wrn = "warning, target "+target+" has equal hits to several queries, cannot decide"
                        warninglist.append(wrn)
                        messagefunc(wrn, debugfile)
    return cond


def bltableout(output, bltableout_file):
    #bltableout_file = open(fname, "a")
    for key, value in sorted(output.items()):
        print >> bltableout_file, key, value
    #bltableout_file.close()

#function to return a list for dict like this:
#dict[query] = [target_f, target_r, target_b, query_f, query_r, query_b, eval, bitscore]
#query - query name, target_f - target start pos, target_r - target end, target_b - forward or reverse target direction
#query_f - query start, query_r - query end, query_b - query direction
def rowfunc(row):
    if int(row[8]) < int(row[9]):
        target_b = True
    else:
        target_b = False
    target_f = int(row[8])
    target_r = int(row[9])
    #check query
    if int(row[6]) < int(row[7]):
        query_b = True
    else:
        query_b = False
    query_f = int(row[6])
    query_r = int(row[7])
    return [target_f, target_r, target_b, query_f, query_r, query_b, float(row[10]), float(row[11])]
    
#function to compute overlap between two ranges supplied as lists with start and end
#returns overlap value
def getOverlap(a, b):
    a0=min(a)
    a1=max(a)
    b0=min(b)
    b1=max(b)
    #print >> debugfile, "interval computation", a, b, a0, a1, b0, b1, max(0, min(a1, b1) - max(a0, b0))
    return max(0, min(a1, b1) - max(a0, b0))

#function for appending the seqeunce to the alignemnt
#returns 1 if success, 0 otherwise
def seqwritefunc(sequence, qname, tname, seqname):
    fhandle = open("./modified/"+qname, "a")
    finalseq = SeqRecord(sequence)
    if contignum == 1 or seqname == "none":
        finalseq.id = tname
    else:
        finalseq.id = tname+"|"+seqname
    finalseq.name =""
    finalseq.description =""
    # if (opt == "-mss" or opt == "-ss") and (float(len(seq.seq)) / loclenfunc(locusfname) > 0.8):
    #     SeqIO.write(seq, fhandle, "fasta")
    #     return 1
    SeqIO.write(finalseq, fhandle, "fasta")
    return 1
    fhandle.close()

# #alternative function for writing, all in 1 file
# def altwritefunc(seq, opt, locusfname, seqname, contignum):
#     fhandle = open("./"+seqname+"_conSeqs.fasta", "a")
#     if contignum == 1:
#         seq.id = locusfname.split("_")[-1].split(".")[0]
#     else:
#         seq.id = locusfname.split("_")[-1].split(".")[0]+"|"+locusfname.split("-copy")[1]
#     seq.name =""
#     seq.description =""
#     print >> fhandle, ">"+seq.id
#     print >> fhandle, seq.seq
#     return 1
#     fhandle.close()

#function for determining the length of the query alignment
#returns length (int)
def loclenfunc(locusfname):
    rhandle = open("./modified/"+locusfname.split("-copy")[0], "r")
    ali = list(SeqIO.parse(rhandle, "fasta"))
    return len(ali[0])
    rhandle.close()

def compute_ranks(hits):
    eval_max = []
    bitscore_avg = []
    coord = []
    for key, hit in hits.items():
        eval_max.append(hit[6])
        bitscore_avg.append(hit[7])
        coord.append(hit[0])
        coord.append(hit[1])
    return [min(eval_max), float(sum(bitscore_avg)) / max(len(bitscore_avg), 1), min(coord), max(coord)]

def hit_sticher(inpdict, extractiontype):
    outlist = []
    #get best item and its direction
    bhit = -1
    best = 1
    for key, val in inpdict.items():
        if bhit == -1:
            bhit = val[7]
        if val[6] < best or (val[6] == best and val[7] > bhit):
            best = val[6]
            bhit = val[7]
            if val[2] == val[5]:
                direct = True
            else:
                direct = False
            if extractiontype == "-n" or extractiontype == "-s" or extractiontype[:2] == "-e" or len(inpdict.keys()) == 1:
                if extractiontype == "-n":
                    outlist = [direct, [-1, -1, val[3], val[4]]]
                else:
                    #Option -e is taken care of at the moment of seq extraction
                    outlist = [direct, [val[0],val[1], val[3], val[4]]]
                    if extractiontype == "-a" or extractiontype == "-b":
                        messagefunc("stiching was not checked, one element", debugfile)
    if extractiontype == "-a" and len(inpdict.keys()) > 1 or extractiontype == "-b" and len(inpdict.keys()) > 1:
        messagefunc("running hit overlapper...", debugfile)
        stichlist = inpdict.values()
        ovlp = True
        while ovlp:
            if len(stichlist) == 1:
                break
            else:
                combos = list(itertools.combinations(range(len(stichlist)), 2))
                messagefunc("iteration length"+str(len(combos))+", "+str(len(stichlist)), debugfile)
                for comb in range(len(combos)):
                    #print >> debugfile, "comparing ...", combos,stichlist
                    tovlp = getOverlap(stichlist[combos[comb][0]][0:2],stichlist[combos[comb][1]][0:2])
                    qovlp = getOverlap(stichlist[combos[comb][0]][3:5],stichlist[combos[comb][1]][3:5])
                    if tovlp > 0 and qovlp > 10: #somehow make the second threshold less arbitrary... perhaps percent?
                        #stichlist[combos[comb][0]] # this is whole record with 8 elements
                        stichlist[combos[comb][0]] = [min(stichlist[combos[comb][0]][0], stichlist[combos[comb][1]][0],stichlist[combos[comb][0]][1], stichlist[combos[comb][1]][1]), max(stichlist[combos[comb][0]][0], stichlist[combos[comb][1]][0],stichlist[combos[comb][0]][1], stichlist[combos[comb][1]][1]), direct, min(stichlist[combos[comb][0]][3], stichlist[combos[comb][1]][3],stichlist[combos[comb][0]][4], stichlist[combos[comb][1]][4]), max(stichlist[combos[comb][0]][3], stichlist[combos[comb][1]][3],stichlist[combos[comb][0]][4], stichlist[combos[comb][1]][4])]
                        del stichlist[combos[comb][1]]
                        messagefunc("overlapped, breaking", debugfile)
                        
                        ovlp = True
                        break
                    elif tovlp == 0 and qovlp > 10: #target non overlapping, but query overlapps - do not stich, remove?
                        #print >> debugfile, "TEST", stichlist[combos[comb][0]],stichlist[combos[comb][1]]
                        if abs(stichlist[combos[comb][0]][0]-stichlist[combos[comb][0]][1]) >= abs(stichlist[combos[comb][1]][0]-stichlist[combos[comb][1]][1]):
                            del stichlist[combos[comb][1]]
                        else:
                            del stichlist[combos[comb][0]]
                        messagefunc("bad hit region, deleting the shortest, breaking", debugfile)
                        ovlp = True
                        break
                    elif tovlp > 10 and qovlp == 0: #target overlapping but query is not - remove as well
                        if abs(stichlist[combos[comb][0]][0]-stichlist[combos[comb][0]][1]) >= abs(stichlist[combos[comb][1]][0]-stichlist[combos[comb][1]][1]):
                            del stichlist[combos[comb][1]]
                        else:
                            del stichlist[combos[comb][0]]
                        messagefunc("bad hit region, deleting the shortest, breaking", debugfile)
                        ovlp = True
                        break
                    else:
                        ovlp = False
                if not ovlp:
                    messagefunc("no more overlaps", debugfile)
        #RUN STICHER and return margins and also sequence of regions and gaps
        messagefunc("running hit sticher...", debugfile)
        median_coords = {}
        start_coords = {}
        end_coords = {}
        start_target = {}
        end_target = {}
        for chunk in range(len(stichlist)):
            #print inplist[chunk][3],inplist[chunk][4]
            median_coords[chunk] = median([stichlist[chunk][3],stichlist[chunk][4]])
            start_coords[chunk] = min(stichlist[chunk][3],stichlist[chunk][4])
            end_coords[chunk] = max(stichlist[chunk][3],stichlist[chunk][4])
            start_target[chunk] = min(stichlist[chunk][0],stichlist[chunk][1])
            end_target[chunk] = max(stichlist[chunk][0],stichlist[chunk][1])
        gapstart = 0
        outlist = [direct]
        for key in sorted(median_coords, key=lambda x: median_coords[x]):
            if gapstart > 0:
                #print >> debugfile, "gap", start_coords[key]-gapstart
                outlist.append(start_coords[key]-gapstart)
            #print >> debugfile, key, median_coords[key], "start", start_coords[key], "end", end_coords[key]
            outlist.append([start_target[key], end_target[key], start_coords[key], end_coords[key]])
            gapstart = end_coords[key]
    return outlist

def median(lst): #taken from https://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def contig_overlap(inplist):
    tab = {}
    for target in inplist:
        #print target[1][-1]
        tab[target[0]] = [min(target[1][1][2],target[1][1][3]),max(target[1][-1][2],target[1][-1][3])]
    flatlist = tab.values()
    messagefunc("checking contig overlap", debugfile)
    
    combos = list(itertools.combinations(range(len(flatlist)), 2))
    ovlp = False
    for comb in range(len(combos)):
        #print combos[comb][1], flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2]
        if getOverlap(flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2]) > 0:
            #messagefunc("overlapping contigs", debugfile)        
            #print "overlapping contigs", flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2]
            ovlp = True
            break
    if ovlp:
        return True
    else:
        return False

def contig_sticher(inplist):
    messagefunc("running contig sticher...", debugfile)
    median_coords = {}
    start_coords = {}
    end_coords = {}
    for target in inplist:
        median_coords[target[0]] = median([min(target[1][1][2],target[1][1][3]),max(target[1][-1][2],target[1][-1][3])])
        start_coords[target[0]] = min(target[1][1][2],target[1][1][3])
        end_coords[target[0]] = max(target[1][-1][2],target[1][-1][3])
    gapstart = 0
    output = []
    for key in sorted(median_coords, key=lambda x: median_coords[x]):
        if gapstart > 0:
            #print >> debugfile, "gap", start_coords[key]-gapstart
            output.append(start_coords[key]-gapstart)
        #print >> debugfile, key, median_coords[key], "start", start_coords[key], "end", end_coords[key]
        output.append(key)
        gapstart = end_coords[key]
    return output

def get_sequence(inplist, seq, extractiontype):
    finalseq = Seq("")
    if extractiontype == "-a":
        for i in inplist[1:]:
            if type(i) is not int:
                start = min(i[0],i[1])
                end = max(i[0],i[1])
                if inplist[0]:
                    finalseq += seq.seq[start:end]
                else:
                    finalseq += seq.seq[start:end].reverse_complement()
            else:
                if i > 0:
                    finalseq += Seq("N"*i)
                else:
                    finalseq += Seq("N")
        #messagefunc("running extractor: final range as is, length "+str(len(finalseq)), debugfile)
    elif extractiontype == "-b":
        f1 = True
        for i in inplist[1:]:
            if type(i) is not int:
                if f1:
                    start = min(i[:2])
                    end = max(i[:2])
                    f1 = False
                else:
                    start = min(start, i[0], i[1])
                    end = max(end, i[0], i[1])
        finalseq = seq.seq[start:end]
        #messagefunc("running extractor: final range "+str([start, end])+", length "+str(len(finalseq)), debugfile)
    elif extractiontype == "-n":
        finalseq = seq.seq
        #messagefunc("running extractor: final range is full, length "+str(len(finalseq)), debugfile)
    elif extractiontype == "-s":
        start = min(inplist[1][0],inplist[1][1])
        end = max(inplist[1][0],inplist[1][1])
        finalseq = seq.seq[start:end]
        #messagefunc("running extractor: final range "+str([start, end])+", length "+str(len(finalseq)), debugfile)
    elif extractiontype[:2] == "-e":
        flank = int(extractiontype[2:])
        start = min(inplist[1][0],inplist[1][1])
        end = max(inplist[1][0],inplist[1][1])
        if start - flank < 0:
            start = 0
        else:
            start = start - flank
        if end + flank > len(seq.seq)-1:
            end = len(seq.seq)-1
        else:
            end = end + flank
        finalseq = seq.seq[start:end]
        #messagefunc("running extractor: final range "+str([start, end])+", length "+str(len(finalseq)), debugfile)
    if not inplist[0] and extractiontype != "-a":
        finalseq = finalseq.reverse_complement()
    return finalseq

def dumper(inplist, extractiontype):
    finalseq = Seq("")
    for i in inplist:
        if type(i) is not int:
            finalseq += i
        else:
            if i > 0 and extractiontype != "-n":
                finalseq += Seq("N"*i)
            else:
                finalseq += Seq("N")
    return finalseq
#---------------------------------------------------------------------

debugfile = open("mainblparser.log", "w")

qout = open("mainblparser_qtable.tab", "w")
tout = open("mainblparser_ttable.tab", "w")

messagefunc("mainblparser run with option "+filefolder+" selected", debugfile, False)
messagefunc("command line parameters: "+' '.join(sys.argv), debugfile, False)


#make modified dir
messagefunc("make modified dir...", debugfile, False)
mkdirfunc()

#copy files
messagefunc("copy files...", debugfile, False)
copyfunc()

#multi db option
if filefolder == "-M":
    #reading the blastfile
    blastlist = glob.glob(blastfilearg+"/*.blast")
    translist = glob.glob(targetf+"/*.fasta")
elif filefolder == "-S":
    blastlist = [blastfilearg]
    translist = [targetf]
else:
    messagefunc("incorrect option @@@"+filefolder, debugfile, False)
    sys.exit()

messagefunc("list of target fasta files detected (mask *.fasta):", debugfile, False)
for l in translist:
    messagefunc(l, debugfile)

#debug vars
number = 0
numberset = set()
totalloci = 0
#parsing blast files
messagefunc("parsing blast files...", debugfile, False)

b1 = 0
for b in blastlist:
    b1 += 1
    messagefunc("target "+str(b1)+" out of "+str(len(blastlist)), debugfile, False)
    output = readblastfilefunc(b, debugfile) #output 0 is query, 1 is target
    final_table = {}
    final_target_table = {}
    for query in output[0].keys():
        #QUERY PROCESSING: first, rank targets by highest eval, also get average bitscore
        messagefunc("Q: "+query, debugfile)
        ranks = [{},{}] #evail is first, bitscore is second
        for target, hits in output[0][query].items():
            ranks_temp = compute_ranks(hits)
            #print >> debugfile, ranks_temp
            #check reciprocy
            if reciprocate == False or reciprocate == True and reciprocator(output[1][target], query, ranks_temp[2], ranks_temp[3], ranks_temp[0],ranks_temp[1], target):
                ranks[0][target] = ranks_temp[0]
                ranks[1][target] = ranks_temp[1]
            else:
                messagefunc("Reciprocator: target "+target+" removed from query "+query,  debugfile)
                
        #print >> debugfile, ranks
        if len(ranks[0]) == 0:
            messagefunc("EMPTY "+query,  debugfile)
            #print >> debugfile, "EMPTY", query
        else:
            sorted_evals = sorted(ranks[0], key=lambda x: ranks[0][x])
            sorted_bits = sorted(ranks[1], key=lambda x: ranks[1][x], reverse=True)
            #GET TARGETS ORDERED AND STICHED
            targets = []
            #print "Q:", query
            for x in range(len(ranks[0])): #using length of ranks, since some contigs are removed due to better hit elsewhere
                if sorted_evals[0] == sorted_bits[0]:
                    messagefunc("best match: "+sorted_evals[0], debugfile)
                    
                else:
                    messagefunc("eval and bit disagree: "+sorted_evals[0]+" and "+sorted_bits[0], debugfile)
                    
                tname1 = sorted_evals.pop(0)
                del sorted_bits[0]
                #SELECT OPTION:
                targets.append([tname1, hit_sticher(output[0][query][tname1], extractiontype)])
            #print >> debugfile, targets
            #CHECK TARGETS FOR OVERLAP
            if interstich:
                if len(targets) > 1:
                    if contig_overlap(targets):
                        messagefunc("contigs overlapping, no contig stiching", debugfile)        
                        stiching_schedule = "none"
                        #CUTTING OFF EXCESS CONTIGS
                        if len(targets) > contignum and contignum > 0:
                            targets = targets[:contignum]
                    else:
                        stiching_schedule = contig_sticher(targets)
                        #ALL will be stiched to just one
                else:
                    messagefunc("only 1 target, no contig stiching", debugfile)
                    stiching_schedule = "none"
            else:
                messagefunc("-IM deactivated, no contig stiching", debugfile)
                stiching_schedule = "none"
                if len(targets) > contignum and contignum > 0:
                    targets = targets[:contignum]

        final_table[query] = [targets, stiching_schedule]
        for t in targets:
            if t[0] in final_target_table:
                final_target_table[t[0]].append(query)
            else:
                final_target_table[t[0]] = [query]
    bltableout(final_table,qout)
    bltableout(final_target_table,tout)

#####-----------------------------------------------------------------------

    messagefunc("scanning the target fasta file...", debugfile, False)
    
    #get the transcriptome filename, matching blast filename
    for t_file in translist:
        if b[:-6].split("/")[-1] in t_file:
            
            inputf = SeqIO.parse(t_file, "fasta")
            
            seqname = b[:-6].split("/")[-1]
            target_db_name = t_file.split("/")[-1]
            break
    #print >> debugfile, "target:", target_db_name, "; target name:", seqname
    if not inputf:
        messagefunc("error, the target fasta file is not found", debugfile, False)
        break
    c1 = len(final_target_table)
    messagefunc("searching for contigs in: "+target_db_name+", total number of contigs: "+str(c1), debugfile, False)

    for seq in inputf: #going over seqs in target file
        if seq.id in final_target_table: #looking up same seq in target file
            for qname in final_target_table[seq.id]: #checking it's queries
                for t in range(len(final_table[qname][0])): #looking for target in the query table
                    if final_table[qname][0][t][0] == seq.id: #found target in the query table
                        if final_table[qname][1] == "none":
                            #extraction
                            messagefunc(str(c1)+" EXTRACTING: contig "+final_table[qname][0][t][0]+", query "+qname, debugfile)
                            s1 = get_sequence(final_table[qname][0][t][1], seq, extractiontype)
                            print >> debugfile, "- EXTRACTING: final seq", s1[:10], "ranges", final_table[qname][0][t][1]
                            seqwritefunc(s1, qname,target_db_name, seq.id)
                        else:
                            s1 = get_sequence(final_table[qname][0][t][1], seq, extractiontype)
                            final_table[qname][1][final_table[qname][1].index(final_table[qname][0][t][0])] = s1
                            messagefunc(str(c1)+" BUCKET: contig "+final_table[qname][0][t][0]+", query "+qname, debugfile)
                            dump_bucket = True
                            for buck1 in final_table[qname][1]:
                                if type(buck1) is str:
                                    dump_bucket = False
                                    break
                            if dump_bucket:
                                s1 = dumper(final_table[qname][1], extractiontype)
                                messagefunc(str(c1)+" EXTRACTING: bucket "+qname+" dumped", debugfile)
                                print >> debugfile, "- EXTRACTING: final seq", s1[:10]#, "ranges", final_table[qname][1]
                                seqwritefunc(s1, qname,target_db_name, "none")
                        #cleanup
                        del final_table[qname][0][t]
                        break # breaking from target table
            #clean up after all qs are done
            del final_target_table[seq.id]
            c1 = c1 - 1
        if len(final_target_table) == 0:
            messagefunc(str(c1)+" search finished", debugfile, False)
            break

print len(warninglist)
for w in warninglist:
    print w

print >> debugfile, "done"
debugfile.close()
qout.close()
tout.close()
print "done"