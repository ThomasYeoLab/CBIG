#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

rep_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}

cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
cmd="${cmd} python ../cbig/Chen2024/CBIG_dnn_transfer_learning.py --src_dataset ${1} --tar_dataset ${2}"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 48:00:00 -mem 256G -ncpus 4 -ngpus 1 \
-name "MMM_transfer_learning" -joberr "${log_dir}/Transfer_learning_err_${1}2${2}.txt" \
-jobout "${log_dir}/Transfer_learning_out_${1}2${2}.txt"
