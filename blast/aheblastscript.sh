echo "erasing the content of the blast output file..."
> blast.blast
########
# argv #
#$1 ahe#
#$2 tri#
#$3 ext#
#$4 clc#
#$5 cex#
#$6 evl#
#$7 num#
########
echo "blasting queries..."
COUNT=0
declare -i CT2=0
declare -i TOTAL1=$(ls $2/*$3 | cat | wc -l) #number of trinity files
declare -i TOTAL2=$(ls $1 | cat | wc -l) #number of AHE files
TOTAL2=$TOTAL2*$TOTAL1*2
for f in $(ls $1)
do
	zo=$(( CT2*100 / TOTAL2 ))
	echo -ne "                                                    \r"
	echo -ne $zo"%\treading... \r"
	while read LINE
	do
		if [[ "{$LINE}" =~ ">" ]]
		then
			echo \>$f > fasextr.blast
			COUNT=1
			continue
		elif [[ $COUNT == 1 ]]
		then
			if [[ "${LINE}" =~ ">" ]]
			then
				COUNT=0
				break
 			else
				echo $LINE >> fasextr.blast
	 		fi
	 	fi
	done < $1/$f
 	for lf in $(ls $2/*$3)
 	do
 		org=$(echo $lf | cut -d"." -f2 | cut -d"/" -f4)
 		echo -ne "                                                                        \r"
 		echo -ne $zo"%\tblast $f against $lf...\r"
 		CT2=$CT2+1
 		tblastx -db $lf -query fasextr.blast -out output.blast -outfmt 6 -evalue $6 -num_threads $7 -num_alignments 1
 		cat output.blast | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge >> blast.blast
 		for clcf in $(ls $4/*$5)
 		do
 			clcorg=$(echo $clcf | cut -d"." -f2 | cut -d"/" -f4 | cut -d"_" -f1)
 			if [[ $clcorg == $org ]]
 			then
 				echo -ne "                                                                          \r"
 				echo -ne $zo"%\tblast $f against $clcf...\r"
 				CT2=$CT2+1
 				tblastx -db $clcf -query fasextr.blast -out output.blast -outfmt 6 -evalue $6 -num_threads $7 -num_alignments 1
 				cat output.blast | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge >> blast.blast
 			fi
 		done
 	done
done
echo -ne "                                                                          \r"
echo "100%..."
echo "done"