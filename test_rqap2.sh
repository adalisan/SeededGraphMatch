# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#

/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "run_sim_experiment_bitflip_rqap2_par_single_job_for_cluster"

echo "" 
echo "Done at " `date` 

