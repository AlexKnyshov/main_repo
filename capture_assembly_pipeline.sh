#$1 - barcode list
#$2 - location of files
module load flash
declare -i ph=1
while read LINE
do
	echo $LINE
	mkdir PH$ph-$LINE
	# pair1=$(find $2 -maxdepth 1 -name "*pair1_$LINE*")
	# pair2=$(find $2 -maxdepth 1 -name "*pair2_$LINE*")
	pair1=$2flowcell692_lane1_pair1_$LINE.fastq.gz
	pair2=$2flowcell692_lane1_pair2_$LINE.fastq.gz
	echo $pair1
	echo $pair2

	cd PH$ph-$LINE
	echo "trimming..."
	java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $pair1 $pair2 output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:64
	echo "read merging..."
	flash output_forward_paired.fq.gz output_reverse_paired.fq.gz
	echo "prep files..."
	gunzip output_forward_unpaired.fq.gz 
	gunzip output_reverse_unpaired.fq.gz 
	cat output_forward_unpaired.fq out.notCombined_1.fastq > forward_up.fq
	cat output_reverse_unpaired.fq out.notCombined_2.fastq > reverse_up.fq
	cat forward_up.fq reverse_up.fq out.extendedFrags.fastq > PH$ph-$LINE.fq
############
	perl ~/tools/prinseq-lite-0.20.4/prinseq-lite.pl -verbose -fastq PH$ph-$LINE.fq -lc_method dust -lc_threshold 45 -out_good PH$ph-$LINE-good
############
	echo "assembly..."
	~/tools/SPAdes-3.10.1-Linux/bin/spades.py -o spades_out-PH$ph -s PH$ph-$LINE-good.fastq

	cd ..
	echo "done"
	ph=$ph+1
done < $1