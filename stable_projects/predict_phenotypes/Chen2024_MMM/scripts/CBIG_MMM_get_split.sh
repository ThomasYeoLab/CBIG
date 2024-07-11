#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.
# rep_dir="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
rep_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}
cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
cmd="${cmd} python ../cbig/Chen2024/CBIG_get_split.py --dataset $1"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 00:30:00 -mem 16G -ncpus 1  -name "MMM_get_split" \
    -joberr "${log_dir}/get_split_${1}_err.txt" -jobout "${log_dir}/get_split_${1}_out.txt"
