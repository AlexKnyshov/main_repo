#PBS -l nodes=1:ppn=16 -l walltime=10:00:00 -l mem=16gb
module load RAxML/8.2.3
cd $PBS_O_WORKDIR
#if [ ! -d "trees" ]; then
#  mkdir trees
#fi
raxmlHPC-PTHREADS-SSE3 -T 16 --silent -f a -c 25 -p 12345 -x 12345 -m GTRCAT -n $(sed -n ${PBS_ARRAYID}p filelist.lst)".tre" -N 100 -s "./phylip/"$(sed -n ${PBS_ARRAYID}p filelist.lst)
