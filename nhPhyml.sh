# $1 - input file
# $2 - tree file
# $3 - number of BS
echo "input file: $1, tree file: $2, number of BS: $3"
~/tools/nhPhyml -sequences=$1 -tree=$2 -format=s -eqfreq=lim -numeqfreq=5
echo "run bootrep.py"
python ~/tools/bootrep.py $1 $3
echo "running bootstrep analyses"
declare -i count=0
for alf in $1.r*
do
	echo "replicate $count, input file $alf"
	~/tools/nhPhyml -sequences=$alf -tree=$2 -format=i -eqfreq=lim -numeqfreq=5
	count=$count+1
done
for tree in $(ls ./*Eq.tree)
do
	echo "converting $tree"
	Rscript ~/tools/tree_conv.R $tree
done
echo "run bootsupport.py"
python ~/tools/bootsupport.py $1_nhPhymlEq
echo "done"