#!/bin/sh
#
# This function performs leave-p-out cross-validation workflow for Elasticnet Regression.
# The function requires 3 arguments.
# 1. param_file: the path to the param file copied and modified from KRR. This is the 
#                output from `CBIG_ICCW_LRR_get_parameters.m`.
# 2. start_fold: an integer of the index of the fold to start the cross validation process from
# 3. fold_number: an integer of the number of folds to process, starting from the start_fold
# 4. LF: path to log file
# 5. root_dir: path to elasticnet scripts
# 
# An example of running this script:
# ./CBIG_ICCW_Elasticnet_job.sh /home/leon_ooi/storage/jianzhong_reliability/results/LRR/behav_1/400 \
#    1 18 /home/leon_ooi/storage/jianzhong_reliability/results/LRR/behav_1/400/behav1_foldStart1.log
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###################################################################################################
# arguments to be read in
###################################################################################################
param_file=$1
start_fold=$2
fold_number=$3
LF=$4
root_dir=$5

###################################################################################################
# run Elasticnet
###################################################################################################
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; load $param_file; \
    CBIG_ICCW_run_Elasticnet_workflow_parallel_folds(param, $start_fold, $fold_number); exit; " >> $LF 2>&1

