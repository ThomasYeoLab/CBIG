#!/bin/sh
#####
# This function is the wrapper script for running regression models that combine multiple features for the
# HCP in our paper.
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

# find main directory
parentdir=$(dirname $(dirname "$(readlink -f "$0")"))
# create ouput folder: MODIFY THIS FOR DIFFERENT OUTPUT DIRECTORY
output_dir=/home/leon_ooi/storage/Multimodal_prediction_project/replication/HCP
if [ ! -d $output_dir ]; then mkdir -p $output_dir; fi

### run second level regressions (stacking and multiKRR)
regression_arr=("stacking" "multiKRR") 
for reg in ${regression_arr[@]}; do
    $parentdir/regression/HCP/utilities/CBIG_MMP_HCP_schedule_regression.sh $reg $output_dir 
done

