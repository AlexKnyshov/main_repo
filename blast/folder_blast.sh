echo "starting..."
if [[ $7 == y ]]
then
	echo "loading modules"
	module load ncbi-blast
else 
	echo "no modules loaded"
fi
echo "erasing the content of the blast output file..."
rm *.blast
########################################
# argv                                 #
#$1 folder with fastas - AHE data      # 
#$2 path to database or folder with dbs#
#$3 threshold                          #
#$4 method: tblastx, blastn            #
#$5 cpu: 4					           #
#$6 list (for multiple dbs)			   #
########################################
if [ -z "$6" ]
  then
    echo "single db blast option selected"
    echo "blasting queries:"
	COUNT=0
	declare -i CT2=0
	declare -i TOTAL=$(ls $1/*.fas | cat | wc -l) #number of ahe files
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
		name=$(echo $f | rev | cut -d'/' -f-1 | rev)
	 	CT2=$CT2+1
	 	echo -ne "                                                                          \r"
	 	echo -ne $zo"% blast $name against $2...\r"
	 	touch output.blast
	 	$4 -db $2 -query fasextr.blast -out output.blast -outfmt 6 -num_threads $5 -evalue $3
	 	blname=$(echo $2 | rev | cut -d'/' -f-1 | rev | cut -d'.' -f1)
	 	cat output.blast | sort -k1,1 -k2,2 -k11,11g -k12,12nr | sort -u -k2,2 --merge | sort -k1,1 -k11,11g -k12,12nr >> $blname.blast
	done
  else
  	echo "multiple db blast option selected, number of dbs: $(cat $6 | wc -l)"
  	echo "blasting queries:"
	COUNT=0
	declare -i CT2=0
	declare -i TOTAL=$(ls $1/*.fas | cat | wc -l)*$(cat $6 | wc -l) #number of files
	zo=$(( CT2*100 / TOTAL ))
	for f in $(ls $1/*.fas)
	do
		
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
		name=$(echo $f | rev | cut -d'/' -f-1 | rev)
		while read fdb
		do
			CT2=$CT2+1
			echo -ne "                                                                          \r"
	 		echo -ne $zo"% blast $name against $fdb...\r"
	 		zo=$(( CT2*100 / TOTAL ))
	 		touch output.blast
	 		$4 -db $2/$fdb -query fasextr.blast -out output.blast -outfmt 6 -num_threads $5 -evalue $3
	 		cat output.blast | sort -k1,1 -k2,2 -k11,11g -k12,12nr | sort -u -k2,2 --merge | sort -k1,1 -k11,11g -k12,12nr >> $fdb.blast
	 	done < $6
	done
fi
echo ""
rm output.blast fasextr.blast
echo "done"