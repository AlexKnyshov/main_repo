echo "starting..."
if [ "$HOSTNAME" = pigeon ]; then
    echo "loading modules"
    module load mafft
elif [ -z "$PBS_ENVIRONMENT" ]; then
	echo "no modules loaded"
else 
	echo "loading modules"
    module load mafft
fi
#$1 - folder with alignments
#$2 - algorithm: linsi, ginsi, einsi
#$3 - direction: adjust, noadjust, slow
#$4 - number of processers
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

for alf in $1*fas
do
	echo "processing file" $alf
	outf=$(echo $alf | cut -d"/" -f3 | cut -d"." -f1,2)
	echo "./realigned/$outf.fas"
	mkdir ./realigned
	mafft $ADJ $ALG --thread $4 --inputorder "$alf" > "./realigned/$outf.fas"
	#mafft --genafpair --maxiterate 1000 --adjustdirection $alf
done