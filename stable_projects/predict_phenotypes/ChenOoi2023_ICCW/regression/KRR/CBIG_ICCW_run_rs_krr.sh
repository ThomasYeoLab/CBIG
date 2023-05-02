#!/bin/bash
# This function runs the kernel regressions in ChenOoi2023.
# When calling this function, sample size should be specified as the first argument.
# Sample size can be 400, 800, 2000, 3000 or 5260.
# The path to the folders containing the results should be specified as the second argument. 
#
# An example of how to call this function as follows:
#     ./CBIG_ICCW_run_rs_krr.sh 400 ~/storage/ICCW_replication
#
# Check the common variables section to see whether the data directories are correct.
# List of directories and files in common variables:
# 1. top_outdir: Output path for the results
# 2. scripts_dir: The folder in which the KRR scripts are stored
# 3. repdata_dir: The folder in which the FC data and covariates are stored
# 4. subject_list: Txt file of subject IDs to run analysis on
# 5. csv_file: csv containing subject IDs, y variables and covariates
# 6. y_list: Txt file of y variables to predict
# 7. covariates: Txt file of covariates to regress from y
# 8. feature_file: Features used for prediction (i.e. FC matrix for all subjects)
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# sample size and results dir to be read in when running script
sample_size=$1
results_dir=$2

###################################################################################################
# set common variables
###################################################################################################
# directories (PLEASE MODIFY)
top_outdir=$results_dir/KRR
scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/KRR
TRBPC_repdata_dir=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
ICCW_repdata_dir=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW
outdir=${top_outdir}/${sample_size}

# files
subject_list=${ICCW_repdata_dir}/input/ICCW_5260_subject_list_unrelated.txt
csv_file=${ICCW_repdata_dir}/input/ICCW_score_motion_mf_w_componentscores.csv
y_list=${ICCW_repdata_dir}/input/ICCW_ABCD_y_variables.txt
covariate_list=${ICCW_repdata_dir}/input/ICCW_covariates_rs.txt
feature_file=${TRBPC_repdata_dir}/FC/FC_subjects_rs_all_score_mf.mat

# settings
FD_file=none
DVARS_file=none
outstem=all_score
curr_task=rs

###################################################################################################
# run single-kernel regression
###################################################################################################
${scripts_dir}/CBIG_ICCW_KRR_LpOCV_workflow.sh -subject_list ${subject_list} \
-feature_file ${feature_file} -outstem ${outstem} -y_list ${y_list} -covariate_list ${covariate_list} \
-FD_file ${FD_file} -DVARS_file ${DVARS_file} -outdir ${outdir} -csv_file ${csv_file} -num_leave_out 5 \
-sample_size ${sample_size}
