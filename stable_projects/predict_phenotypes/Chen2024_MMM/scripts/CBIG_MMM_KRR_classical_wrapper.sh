#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

# Initialize directory
input_dir="$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Chen2024_MMM"
base_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM"
code_dir="${base_dir}/KRR_CLASSICAL"
rep_dir="${base_dir}/replication"
dataset=$1
# Initialize input and output
ks=(10 20 50 100 200)
final_phe_list="${input_dir}/${1}/${1}_phe_list.txt"
phe_csv="${input_dir}/${1}/${1}_phe_tab.csv"
subj_list="${input_dir}/${1}/${1}_subj_list.txt"
fc_mat="${input_dir}/${1}/${1}_fc.mat"
output="${rep_dir}/output_KRR_classical_${1}"

# Run KRR by submitting job
for k in "${ks[@]}"; do
  while read p; do
    n_job="$(qstat -u the | wc -l)"
    while [ $n_job -gt 50 ]; do
      echo $n_job large than 50, sleep for 1 mins
      sleep 1m
      n_job="$(qstat -u the | wc -l)"
    done
    sh ${base_dir}/scripts/CBIG_MMM_KRR_classical_submit_job.sh $p $output $phe_csv $code_dir $k \
    $subj_list $fc_mat $dataset
  done <$final_phe_list
done
