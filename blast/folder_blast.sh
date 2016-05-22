echo "erasing the content of the blast output file..."
> blast.blast
########
# argv #
#$1 folder with fastas# 
#$2 path to query#
#$3 threshold
########
echo "blasting queries..."
COUNT=0
declare -i CT2=0
declare -i TOTAL=$(ls $1/*.fas | cat | wc -l) #number of trinity files
for f in $(ls $1/*.fas)
do
	zo=$(( CT2*100 / TOTAL ))
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
	done < $f
 	echo -ne $zo"% blast $f against $2...\r"
 	CT2=$CT2+1
 	touch output.blast
 	tblastx -db $2 -query fasextr.blast -out output.blast -outfmt 6 -num_threads 4 -num_alignments 1 -evalue $3
 	cat output.blast | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge >> blast.blast
done
echo -ne "                                                                          \r"
echo "100%..."
echo "done"