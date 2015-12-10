#this file describes the work


#starting folders - only ./fasta
#all other folders are being created when scripts are executed

#1) preselect taxa; output is written to ./reduced
python TaxonParser

#2) realign fasta files; output is written to ./realigned
python realignment.py

#3) convert fasta to phylip; output is written to ./phylip
python fstphil.py ./realigned fasta realigned.fas ./phylip phylip-relaxed phylip

#4) compute distances; in last script version -gtr option stopped working on the cluster, but not enough time to tweak it. Should work on any laptop with raxml installed. -p and -115 should work everywhere. Computation takes no more than 15 minutes
python distnew -gtr ./phylip

#5) filter based on distances; output >> filteredloci.tab
python treefilter.py locigtrn.tab
#for min and max we indicated 7.0 and 41.5

#6) run a raxml on the cluster through qsub arrayjob; works only on cluster
#filteredloci.tab should be present in the folder
#instead of 390 user puts how many loci left after filtering (wc -l filteredloci.tab)
qsub -t 1-390 treegen.sh

#7) move trees
bash treselect.sh

#8) process trees
python treeprocsep.py -75 ./trees

#9) filter based on trees; output >> filteredloci.tab
python treefilter.py treesfilter.tab
#for min we indicated 0.35, for max - 0.55

#10) concatenate the remaining loci; output is COMBINED.nex
python concatenator.py ./phylip ./nexus filteredloci.tab

#11a) make COMBINED.phy and partitionfinder config
python partextractor.py

#11b) make partition file for raxml (partitions.prt)
#conversion of COMBINED.nex to phy should be done by other program, i.e. Mesquite
python raxmlpartextractor.py

#12) run the final analysis
raxmlHPC-PTHREADS-SSE3 -T 4 -f a -c 25 -p 12345 -x 12345 -m GTRCAT -n TEST -q partitions.prt -s COMBINED.phy -N 1000
#can be submitted on a cluster as well, but no array job is needed since only 1 dataset is present

#13) render the final tree in FigTree