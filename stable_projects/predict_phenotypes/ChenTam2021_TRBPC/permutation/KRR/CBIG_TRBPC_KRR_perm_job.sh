#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
model=$1
input_dir=$2
outstem=$3
score=$4
perm_start=$5
N_per_job=$6
group=$7
outdir=$8
scripts_dir=`dirname "$(readlink -f "$0")"`

LF=${outdir}/logs/perm_test_score${score}_perm_start${perm_start}.txt
date >> $LF

matlab -nodesktop -nosplash -nodisplay -r " addpath $scripts_dir; \
CBIG_TRBPC_compute_${model}_perm_stats('${input_dir}','${outstem}','${score}',${perm_start},\
${N_per_job},'${group}','${outdir}'); exit; " >> $LF
