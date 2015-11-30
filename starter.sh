qsub -l nodes=1:ppn=8
python treeanalysis.py ./original_data 50
#echo ‘myscript.sh prot.faa 2’ | qsub