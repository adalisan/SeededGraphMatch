# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#
MATLABPATH=/usr/local/bin
FILEPATH=/cis/home/sadali/my_documents/ 
M_FILENAME=run_sim_experiment_bitflip_par_single_job_for_cluster.m
# Name the job #$ -N matlabScript #
cd $FILEPATH 
$MATLABPATH/matlab -nosplash -nodisplay -singleCompThread <<EOF
$M_FILENAME
EOF
echo "" 
echo "Done at " `date`
