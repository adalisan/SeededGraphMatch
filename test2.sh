# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#

/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "run_sim_experiment_bitflip_par_single_job_for_cluster_test"
echo "" 
echo "Done at " `date` 

