########
AHE + transcriptome pipeline:

1) prepare folder with AHE alignments in fasta format
	- to convert into fasta from phylip:
	python fconv.py -a ./PHYLIP_FOLDER phylip-relaxed .phylip ./FASTA_FOLDER fasta .fas
	- to subset taxa in AHE alignments:
	create a list with taxon names
		- to leave specified taxa:
		python removeTaxa.py ./FASTA_FOLDER -a LIST
		- to exclude specified taxa
		python removeTaxa.py ./FASTA_FOLDER -e LIST
		- to leave and rename specified taxa:
		python removeTaxa.py ./FASTA_FOLDER -a LIST #list must contain old and new names separated with comma (oldname,newname)
	- to make the top sequence the longest one (which is good for blast searches):
	python taxon_regroup.py -seqlen ./FASTA_FOLDER
2) prepare transcriptomes - they need to be in one folder, with databases created in the same folder
3) run blast searches
	- against a single database:
	bash folder_blast.sh ./AHE_FOLDER ./TRANSCRIPTOME_FILENAME threshold METHOD[i.e., tblastx] CPU[i.e. 4]
	- against a group of databases:
	bash folder_blast.sh ./AHE_FOLDER ./TRANSCRIPTOME_FOLDER threshold METHOD[i.e., tblastx] CPU[i.e. 4] LIST
	list must specify transcriptome filenames in the folder indicated. if all files in a certain folder shoudl be used,
	this list can be obtained by ls * > list.txt
4) extract found contigs:
	- to extract from a single transcriptome:
	python alexblparser.py BLAST_OUTPUT_FILE TRANSCRIPTOME_FILENAME ./AHE_FOLDER threshold OPTION
	- to extract from a group of transcriptomes:
	python alexblparser.py BLAST_OUTPUT_FOLDER TRANSCRIPTOME_FOLDER ./AHE_FOLDER threshold OPTION
	note that the script will search for all *.blast files in the indicated blast folder
	
	loci are exported to "modified" folder

	OPTIONs available are:
	-n (extracts an entire contig, single transcriptome mode)
	-s (extracts the hit region of a contig, single transcriptome mode)
	-ss (extracts the hit region of a contig, hits shorter than 80% length of query are discarded, single transcriptome mode)
	-mn (extracts an entire contig, group of transcriptomes mode)
	-ms (extracts the hit region of a contig, group of transcriptomes mode)
	-mss (extracts the hit region of a contig, hits shorter than 80% length of query are discarded, group of transcriptomes mode)
	-me (extracts the region of a conting that matches the length of the query, group of transcriptomes mode)
5) realign the loci
	bash align.sh ./FASTA_FOLDER ALGORITHM DIRECTION_OPTION CPU(number)
	algorithm options are: ginsi einsi linsi
	direction options are: adjust, noadjust, slow
	files output to "realigned" folder
6) if adjustment of sequence direction was performed (options adjust or slow), remove R_ prefixes:
	python removeR.py ./FASTA_FOLDER -m
	existing files are modified
7) if trimming is necessary
	python customtrim.py ./FASTA_FOLDER OPTION
	files output to "trimmed"
	available options are:
	-a (trim external gaps until all AHE loci have no missing data)
	-1 (trim external gaps until one AHE locus has no missing data)
	-% (trim external gaps until 90% of sequences have no missing data)
	-refine (trim external gaps until all sequences have no missing data for more than 20 bases)
	-d (trim all positions with missing data)
8) if subsetting based on taxa presence is necessary
	python AHE_taxatrim.py ./FASTA_FOLDER LIST
	a locus will be removed if no listed taxa are present
	files ouput to "subset" folder
