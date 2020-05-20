mkdir upho_trees
for f in $1/RAxML_bipartitions.*
do
	bname=$(echo $f | rev | cut -f1,2 -d. | rev)
	cp $f upho_trees/$bname
done