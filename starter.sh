#qsub -l nodes=1:ppn=8
qsub -S $(which python) -l nodes=1:ppn=8 treeanalysis.py ./original_data 50
#echo ‘myscript.sh prot.faa 2’ | qsub