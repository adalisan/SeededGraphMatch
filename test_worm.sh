# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 



/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "addpath('~/projects/SeededGraphMatch/src/');addpath('~/projects/SeededGraphMatch/lib/');run_worm_experiment_par_single_job_for_cluster"

echo "" 
echo "Done at " `date` 

