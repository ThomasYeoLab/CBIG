#!/bin/bash
#
# This function runs the linear ridge regressions in ChenOoi2023.
# When calling this function, sample size should be specified as the first argument. 
# Sample size can be 400, 800, 2000, 3000 or 5260.
# The path to the folders containing the results should be specified as the second argument. 
#
# An example of how to call this function as follows:
#     ./CBIG_ICCW_LRR_wrapper.sh 400 ~/storage/ICCW_replication
#
# Check the common variables section to see whether the data directories are correct.
# List of directories and files in common variables:
# 1. root_dir: The folder in which the LRR scripts are stored
# 2. util_dir: The folder in which the Elasticnet scripts are stored
# 3. KRR_dir: The folder in which the KRR results are stored
# 4. out_dir: Output path for the results
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# sample size variable to be read in when running script
sample_size=$1
results_dir=$2

###################################################################################################
# set common variables
###################################################################################################
# directories (PLEASE MODIFY)
root_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/LRR
util_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/utilities
KRR_dir=$results_dir/KRR/$sample_size
out_dir=$results_dir/LRR/$sample_size

###################################################################################################
# run LRR
###################################################################################################
mkdir -p $out_dir/logs
for behav in `seq 1 1` # only run for 1 behavior, for full set replace with `seq 1 39`
do 
    # prepare parameters for LRR
    LF=$out_dir/logs/Behav${behav}.log
    matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; CBIG_ICCW_LRR_get_parameters( \
       '$KRR_dir', '$out_dir', $behav); exit; " >> $LF 2>&1
    param_file=$out_dir/behav_${behav}/param.mat

    # run LRR in batches of 18 folds each
    for fold_start in `seq 1 18 252`
    do
        logfile=${out_dir}/logs/behav${behav}_foldStart${fold_start}.log
        cmd="${util_dir}/CBIG_ICCW_Elasticnet_job.sh $param_file $fold_start 18 $logfile $util_dir"
    ${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 24:00:0 -mem 24gb -joberr ${logfile} -jobout ${logfile}\
     -cmd "${cmd}" -name ABCD_LRR
    done
done
