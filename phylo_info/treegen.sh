#PBS -l nodes=1:ppn=16 -l mem=16gb
module load RAxML/8.2.3
cd $PBS_O_WORKDIR
raxmlHPC-PTHREADS-SSE3 -T 16 --silent -f a -c 25 -p 12345 -x 12345 -m GTRCAT -n $(sed -n ${PBS_ARRAYID}p filteredloci.tab)".tre" -N 1000 -s "./phylip/"$(sed -n ${PBS_ARRAYID}p filteredloci.tab)
