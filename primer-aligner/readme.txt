Primer-aligner manual:

- subject seqeunce need to be in fasta format, single seqeunce tested only. Umbiguous bases are not supported.
- primers need to be in fasta format in one file. Reverse primerse will be reversed and complemented as the script runs
- perform blast search as follows:
  blastn -query primers.fas -subject subject.fas -outfmt 6 -task blastn-short | sort -k1,1 -k2,2 -k11,11g -k12,12nr | sort -u -k1,1 > blast.blast
- run the python script as follows:
  python primer-aligner.py subject.fas primers.fas blast.blast -c
