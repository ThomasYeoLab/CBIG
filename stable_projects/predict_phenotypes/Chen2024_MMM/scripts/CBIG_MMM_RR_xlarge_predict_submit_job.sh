#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

rep_dir="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}
cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
cmd="${cmd} python ../cbig/Chen2024/CBIG_rr_xlarge_predict.py --src_dataset $1"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 6:00:00 -mem 128G -ncpus 8 -name "MMM_RR_pred_XL" \
-joberr "${log_dir}/$1_RR_predict_err.txt" -jobout "${log_dir}/$1_RR_predict_out.txt"
