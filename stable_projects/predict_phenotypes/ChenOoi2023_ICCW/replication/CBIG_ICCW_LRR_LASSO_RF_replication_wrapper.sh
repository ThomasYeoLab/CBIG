#!/bin/bash

# This script runs the prediction replication for the 400 and 800 sample sizes (as running for all 
# sample sizes will take too long). The script only runs the LRR, LASSO and RF.
#
# Example:
#     sh CBIG_ICCW_LRR_LASSO_RF_replication_wrapper.sh ~/storage/CBIG_ICCW_replication
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
results_dir=$1
if [ ! -d $results_dir ]; then
	mkdir -p $results_dir
else
	echo "Directory for results_dir already exists!"
fi
scripts_dir=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression
sample_sizes=(400 800) # only run for 400 and 800, if not it will take too long

####################################################
# Submit jobs for LRR
####################################################
for sample in ${sample_sizes[@]}; do
	sh $scripts_dir/LRR/CBIG_ICCW_LRR_wrapper.sh $sample $results_dir
done

####################################################
# Submit jobs for LASSO
####################################################
for sample in ${sample_sizes[@]}; do
	sh $scripts_dir/LASSO/CBIG_ICCW_LASSO_wrapper.sh $sample $results_dir
done

####################################################
# Submit jobs for RF
####################################################
for sample in ${sample_sizes[@]}; do
	sh $scripts_dir/RF/CBIG_ICCW_submit_RF_train.sh $sample $results_dir
done