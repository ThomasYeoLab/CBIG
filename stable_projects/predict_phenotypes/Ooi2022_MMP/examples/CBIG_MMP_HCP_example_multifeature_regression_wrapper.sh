#!/bin/sh
#####
# This function is the wrapper script for running regression models that combine multiple features for the
# HCP in our paper. Please specify the output directory as the first argument.
# Example:
# ./CBIG_MMP_HCP_example_multifeature_regression_wrapper.sh $output_dir
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

# find main directory
parentdir=$(dirname $(dirname "$(readlink -f "$0")"))
output_dir=$1
if [ ! -d $output_dir ]; then mkdir -p $output_dir; fi

### run second level regressions (stacking and multiKRR)
regression_arr=("stacking" "multiKRR") 
for reg in ${regression_arr[@]}; do
    $parentdir/examples/scripts/utilities/CBIG_MMP_HCP_schedule_regression_example.sh \
    $reg $output_dir 
done

