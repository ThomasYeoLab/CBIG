#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################
input_dir=$1
feature_file=$2
sub_fold=$3
score_ind=$4
outdir=$5
model=$6

########## compute the PFM
LF="$outdir/logs/compute_${model}_PFM.log"
date >> $LF
scripts_dir=`dirname "$(readlink -f "$0")"`
matlab -nodesktop -nosplash -nodisplay -r " addpath $scripts_dir; CBIG_TRBPC_compute_${model}_PFM( \
'$feature_file', '$input_dir', '$sub_fold', '$score_ind', '$outdir' );\
exit; " >> ${outdir}/logs/${model}_PFM_score_${score_ind}.txt