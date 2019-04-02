echo "starting..."
if [[ $5 == y ]]
then
	echo "loading modules"
	module load mafft
else 
	echo "no modules loaded"
fi
#$1 - folder with alignments
#$2 - algorithm: linsi, ginsi, einsi
#$3 - direction: adjust, noadjust, slow
#$4 - number of processers
#$5 -y/n, loading modules
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
ADJ=""
if [[ $3 == adjust ]]
then
	ADJ="--adjustdirection"
elif [[ $3 == noadjust ]]
then
	ADJ=""
elif [[ $3 == slow ]]
then
	ADJ="--adjustdirectionaccurately"
fi
mkdir ./realigned/
for alf in $1*fas
do
	echo "processing file" $alf
	outf=$(echo $alf | rev | cut -d"/" -f1 | rev)
	echo "./realigned/$outf"
	mafft $ADJ $ALG --thread $4 --inputorder "$alf" > "./realigned/$outf"
	#mafft --genafpair --maxiterate 1000 --adjustdirection $alf
done