#arguments
#$1 - folder with trinity
#$2 - extension of trinity files (with dot)
#$3 - folder with clc
#$4 - extension of clc files (with dot)

#trinity dbs
for qf in $(ls $1/*$2)
do
	echo $qf
	makeblastdb -in $qf -dbtype nucl -parse_seqids
done
#clc dbs
for clf in $(ls $3/*$4)
do
	echo $clf
	makeblastdb -in $clf -dbtype nucl -parse_seqids
done
