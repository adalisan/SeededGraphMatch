# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#

/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "addpath(strcat(pwd,'/src'));addpath(strcat(pwd,'/lib'));run_sim_experiment_bitflip_hybrid_par_single_job_for_cluster"

echo "" 
echo "Done at " `date` 

