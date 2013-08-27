# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 


/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread  -r "addpath('/cis/home/sadali/my_documents/src/');addpath('/cis/home/sadali/my_documents/lib/');run_enron_experiment_par_single_job_for_cluster"

echo "" 
echo "Done at " `date` 

