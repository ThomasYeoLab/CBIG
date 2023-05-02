#!/bin/bash

# This script runs the predictive feature interpretation jobs that require job submission, and calculates
# the t-statistics, which is consistent across all models. Original regression weights are calulcated for
# KRR and conditional variable importance is calculated for RF.
#
# Example:
#     sh CBIG_ICCW_interpretation_replication_wrapper.sh ~/storage/ICCW_replication
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
# Submit jobs for KRR
####################################################
for sample in ${sample_sizes[@]}; do
	sh $scripts_dir/KRR/CBIG_ICCW_run_krr_weights_job.sh $sample $results_dir
done

####################################################
# Submit jobs for RF
####################################################
for sample in ${sample_sizes[@]}; do
	sh $scripts_dir/RF/CBIG_ICCW_submit_RF_cperm.sh $sample $results_dir
done

####################################################
# Run matlab script for tstats
####################################################
mkdir -p $results_dir/analysis/logs
LF="$results_dir/analysis/logs/tstats.log"
if [ -f $LF ]; then rm $LF; fi
root_dir=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/analysis/all_models
matlab -nodesktop -nosplash -nodisplay -r " addpath ${root_dir}; \
CBIG_ICCW_compute_tstats('$results_dir/KRR', '$results_dir/analysis'); \
exit; " >> $LF 2>&1
