# 
#$ -wd  /cis/home/sadali/my_documents/logs/
#$ -j y 
#$ -S /bin/bash 
#$
#

/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "addpath('/cis/home/sadali/my_documents/src/');addpath('/cis/home/sadali/my_documents/lib/');run_sim_experiment_bitflip_par_single_job_for_cluster"

echo "" 
echo "Done at " `date` 

