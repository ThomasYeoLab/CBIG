#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.
# rep_dir="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
rep_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}
cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
cmd="${cmd} python ../cbig/Chen2024/CBIG_mm_stacking.py --log_stem $1 --tar_dataset $2"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 3:00:00 -mem 64G -ncpus 8  -name "MM" \
    -joberr "${log_dir}/${1}_${2}_err.txt" -jobout "${log_dir}/${1}_${2}_out.txt"
