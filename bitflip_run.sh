# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#
PATH_TO_B=bitflip_sim.sh
echo Calling b.sh

for i in `seq 1 500`; do
qsub $PATH_TO_B
done
echo "" 
echo "Done at " `date` 

