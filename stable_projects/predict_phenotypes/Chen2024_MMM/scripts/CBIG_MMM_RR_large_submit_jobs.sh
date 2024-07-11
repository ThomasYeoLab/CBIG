#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

rep_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}
rep_data_dir="$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/"
src_dataset=$1
layer=$2
phe_list_path="${rep_data_dir}/${src_dataset}/${src_dataset}_phe_list.txt"
num_phe=$(wc -l < ${phe_list_path})

for phe in $(seq 0 $((num_phe - 1))); do
    cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
    cmd="${cmd} python ../cbig/Chen2024/CBIG_rr_large.py --phe_idx ${phe} \
        --src_dataset ${src_dataset} --${layer}"
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 6:00:00 -mem 32G -ncpus 4 -name "MMM_RR_L" \
        -joberr "${log_dir}/${1}_LRR_${2}_${phe}_err.txt" -jobout "${log_dir}/${1}_LRR_${2}_${phe}_out.txt"
done
