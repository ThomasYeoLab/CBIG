#!/bin/bash
# This shell script calls the matlab script to compute regression weights for KRR.
#
# The following inputs are required..
# List of directories and files in common variables:
# 1. krr_input: Features used for prediction (i.e. FC matrix for all subjects)
# 2. krr_result_dir: Directory in which KRR results were saved
# 3. score_ind: A scalar. Index of behavior score to extract weights for
# 4. outdir: output directory for results to be saved in
# 5. scripts_dir: The folder in which the KRR scripts are stored
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###################################################################################################
# read inputs
###################################################################################################
krr_input=$1
krr_result_dir=$2
score_ind=$3
outdir=$4
###################################################################################################
# call matlab function
###################################################################################################
root_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/KRR
LF=$outdir/logs/krr_weights_${score_ind}.log
mkdir -p $outdir/logs
matlab -nodesktop -nosplash -nodisplay -r " addpath $root_dir; \
CBIG_ICCW_compute_krr_weights('$krr_input','$krr_result_dir',$score_ind,'$outdir'); exit; " >> $LF 2>&1
