#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

FC_mean_file=$1
y_file=$2
krr_folds=$3
perm_seed_start=$4
perm_num=$5
outdir=$6
group=$7

scripts_dir=`dirname "$(readlink -f "$0")"`

LF=${outdir}/logs/PFM_perm_start_${perm_seed_start}.txt
date >> $LF
matlab -nodesktop -nosplash -nodisplay -r " addpath $scripts_dir; CBIG_TRBPC_compute_PFM_permutation_stats( \
'$FC_mean_file','$y_file','$krr_folds','$perm_seed_start','$perm_num','$group','$outdir'); \
exit; " >> $LF
