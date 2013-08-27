# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#

/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "addpath(strcat(pwd,'/src'));run_sim_experiment_bitflip_hybrid_par_single_job_for_cluster_t"

echo "" 
echo "Done at " `date` 

