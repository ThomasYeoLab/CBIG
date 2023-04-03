#!/bin/bash

# Example:
#     sh CBIG_GradPar_replication_wrapper.sh ~/storage/Temporary/CBIG_GradPar_replication
#
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
out_dir=$1
echo $out_dir
mkdir -p $out_dir
out_dir=`realpath $1`
out_LRR_dir=$out_dir/LRR
out_KRR_dir=$out_dir/KRR
out_opt_res_dir=$out_dir/opt_res
scripts_dir=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/replication
####################################################
# Submit jobs for KRR
####################################################
sh $scripts_dir/CBIG_GradPar_KRR_each_resolution_submit.sh $out_KRR_dir

####################################################
# Submit jobs for LRR
####################################################
sh $scripts_dir/CBIG_GradPar_LRR_frac_each_resolution_submit.sh $out_LRR_dir

####################################################
# Submit jobs for optimizing resolution parameters
####################################################
# To optimize the resolution, we require the user to generate the prediction results for different
# resolutions first. Here we use the prediction results provided in the reference folder to save time.

mkdir -p $out_opt_res_dir/logs
log_dir=$out_opt_res_dir/logs/opt_res.log

ref_LRR_dir=$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/ref_results/LRR
ref_KRR_dir=$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/ref_results/KRR

cmd="cd ${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Kong2023_GradPar/replication;"
cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
cmd="$cmd \"clear;clc;close all;CBIG_GradPar_optimize_res_wrapper \
${ref_LRR_dir} ${ref_KRR_dir} ${out_opt_res_dir}; exit; \" "
cmd="$cmd | tee -a $log_dir"

ERR_FILE_PATH="${out_opt_res_dir}/job_err_out/job.err"
OUT_FILE_PATH="${out_opt_res_dir}/job_err_out/job.out"
mkdir -p ${out_opt_res_dir}/job_err_out

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 3:00:00 -mem 16G \
-name "GradPar_opt_res" -joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

cmd=""