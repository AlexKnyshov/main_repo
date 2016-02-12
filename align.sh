echo "starting..."
#module load mafft
#$1 - folder with alignments
for alf in $1*fas
do
	echo "processing file" $alf
	outf=$(echo $alf | cut -d"/" -f3 | cut -d"." -f1,2)
	echo "./realigned/$outf.fas"
	mkdir ./realigned
	mafft --adjustdirection --genafpair  --maxiterate 1000 --inputorder "$alf" > "./realigned/$outf.fas"
	#mafft --genafpair --maxiterate 1000 --adjustdirection $alf
done