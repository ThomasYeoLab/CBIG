#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Initialize working directory
input_dir="$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/He2022_MM"
base_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM"
code_dir="${base_dir}/KRR_MM"
rep_dir="${base_dir}/replication"

# Run KRR for phe except sex
split="test"
rngs=1
phe_list="${input_dir}/ukbb_train_final_phe_list.txt"
output="${rep_dir}/output_KRR_mm"
phe_csv="${input_dir}/ukbb_train_${split}_final.csv"
subj_list="${input_dir}/ukbb_train_${split}_subj_list.txt"
fc_mat="${input_dir}/ukbb_train_${split}_pfc.mat"
subj_list_train="${input_dir}/ukbb_train_subj_list.txt"
subj_list_extra="${input_dir}/ukbb_${split}_subj_list.txt"

for j in `seq 1 $rngs`; do
    n_job="$(qstat -u the | wc -l)"
    while [ $n_job -gt 100 ]; do
        echo $n_job large than 100, sleep for 1 mins
        sleep 1m
        n_job="$(qstat -u the | wc -l)"
    done
    sh CBIG_MM_KRR_MM_submit_job.sh $output $phe_csv $code_dir $split $j \
        $phe_list $subj_list $fc_mat $subj_list_train $subj_list_extra
done

# Run KRR for sex
phe_list="31-0.0"
output="${rep_dir}/output_KRR_mm_sex"

for j in `seq 1 $rngs`; do
    n_job="$(qstat -u the | wc -l)"
    while [ $n_job -gt 100 ]; do
        echo $n_job large than 100, sleep for 1 mins
        sleep 1m
        n_job="$(qstat -u the | wc -l)"
    done
    sh CBIG_MM_KRR_MM_submit_job.sh $output $phe_csv $code_dir $split $j \
        $phe_list $subj_list $fc_mat $subj_list_train $subj_list_extra
done