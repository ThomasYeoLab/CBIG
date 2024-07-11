#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

rep_dir="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}
src_dataset=$1
rep_data_dir="/mnt/isilon/CSC2/Yeolab/Data/UKBioBank/CodeMaintanance/ReplicationData/"
phe_list_path="${rep_data_dir}/stable_projects/predict_phenotypes/Chen2024_MMM/${src_dataset}_phe_list.txt"
num_phe=$(wc -l < ${phe_list_path})

for phe in $(seq 0 $((num_phe - 1))); do
    cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
    cmd="${cmd} python ../cbig/Chen2024/CBIG_rr_xlarge_train.py --phe_idx ${phe} --src_dataset ${src_dataset}"
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 72:00:00 -mem 128G -ncpus 8 -name "MMM_RR_train_XL" \
        -joberr "${log_dir}/${1}_RR_${phe}_err.txt" -jobout "${log_dir}/${1}_RR_${phe}_out.txt"
done
