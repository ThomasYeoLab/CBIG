#!/bin/sh
#####
# This script calls the matlab function to run the Haufe inversion for KRR models. 
# User needs to provide the following variables.
# 1. input_dir: The directory in which the brain imaging features are saved.
# 2. results_dir: The directory in which the regression results are results are saved.
# 3. feature:  The outstem of the model to invert (e.g. features_rs).
# 
# EXAMPLE: 
#    CBIG_MMP_HCP_Haufe.sh $input_dir $results_dir $feature 
# EXAMPLE OF HOW TO CALL FUNCTION:
#     CBIG_MMP_HCP_Haufe.sh data_dir/features output_dir 'features_rs' 
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# set directories
script_dir=$(dirname "$(readlink -f "$0")")

# set params
input_dir=$1
results_dir=$2
feature=$3

# Create log file and save params
output_dir=$results_dir/interpretation/$feature
if [ ! -d $output_dir ];then mkdir -p $output_dir; fi
LF="$output_dir/${feature}_interpretation.log"
if [ -f $LF ]; then rm $LF; fi
echo "input_dir = $input_dir" >> $LF
echo "results_dir = $results_dir" >> $LF
echo "feature = $feature" >> $LF

# Call matlab function
matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_HCP_Haufe( \
   '$input_dir','$results_dir', '$feature'); catch ME; display(ME.message); \
   end; exit; " >> $LF 2>&1