echo "starting..."
#module load mafft
#$1 - folder with alignments
#$2 - algorithm: linsi, ginsi, einsi
ALG=""
if [[ $2 == einsi ]]
then
	ALG="--genafpair  --maxiterate 1000"
elif [[ $2 == ginsi ]]
then
	ALG="--globalpair --maxiterate 1000"
elif [[ $2 == linsi ]]
then
	ALG="--localpair --maxiterate 1000"
fi
#echo $ALG

for alf in $1*fas
do
	echo "processing file" $alf
	outf=$(echo $alf | cut -d"/" -f3 | cut -d"." -f1,2)
	echo "./realigned/$outf.fas"
	mkdir ./realigned
	mafft --adjustdirection $ALG --inputorder "$alf" > "./realigned/$outf.fas"
	#mafft --genafpair --maxiterate 1000 --adjustdirection $alf
done