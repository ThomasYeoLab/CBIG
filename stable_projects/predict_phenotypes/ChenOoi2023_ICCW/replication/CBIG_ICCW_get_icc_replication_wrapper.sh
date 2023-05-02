#!/bin/bash

# This script runs completes calculating the predictive feature importance for the Haufe transformed weights
# and calculates the agreement and final prediction accuracy from all models.
#
# Example:
#     sh CBIG_ICCW_get_icc_replication_wrapper.sh ~/storage/ICCW_replication
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
results_dir=$1
if [ ! -d $results_dir ]; then
	echo "results_dir does not exist!"
	exit
fi
scripts_dir=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/analysis
output_dir=$results_dir/analysis
mkdir -p $results_dir/analysis/logs

####################################################
# Get Haufe and ICC for KRR
####################################################
LF="$results_dir/analysis/logs/KRR_acc_icc.log"
if [ -f $LF ]; then rm $LF; fi
matlab -nodesktop -nosplash -nodisplay -r " addpath ${scripts_dir}/KRR; \
CBIG_ICCW_compute_KRR_acc_icc('$results_dir/KRR', 'weights', '$results_dir/analysis'); \
exit; " >> $LF 2>&1

####################################################
# Get Haufe and ICC for LRR
####################################################
LF="$results_dir/analysis/logs/LRR_acc_icc.log"
if [ -f $LF ]; then rm $LF; fi
matlab -nodesktop -nosplash -nodisplay -r " addpath ${scripts_dir}/LRR_LASSO; \
CBIG_ICCW_compute_Elasticnet_acc_icc('LRR', '$results_dir', '$results_dir/analysis'); \
exit; " >> $LF 2>&1

####################################################
# Get Haufe and ICC for LASSO
####################################################
LF="$results_dir/analysis/logs/LASSO_acc_icc.log"
if [ -f $LF ]; then rm $LF; fi
matlab -nodesktop -nosplash -nodisplay -r " addpath ${scripts_dir}/LRR_LASSO; \
CBIG_ICCW_compute_Elasticnet_acc_icc('LASSO', '$results_dir', '$results_dir/analysis'); \
exit; " >> $LF 2>&1

####################################################
# Get Haufe and ICC for RF
####################################################
LF="$results_dir/analysis/logs/RF_acc_icc.log"
if [ -f $LF ]; then rm $LF; fi
matlab -nodesktop -nosplash -nodisplay -r " addpath ${scripts_dir}/RF; \
CBIG_ICCW_compute_RF_acc_icc('$results_dir', '$results_dir/analysis'); \
exit; " >> $LF 2>&1