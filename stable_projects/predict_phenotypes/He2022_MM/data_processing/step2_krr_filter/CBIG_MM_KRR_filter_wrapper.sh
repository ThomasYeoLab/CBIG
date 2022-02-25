#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ARG1=${1:-normal}

# Initialize working directory
USERNAME=`whoami`
base_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM"
stage_1_out_dir="${base_dir}/data_processing/step1_coarse_filter/output"
stage_2_dir="${base_dir}/data_processing/step2_krr_filter"
stage_3_out_dir="${base_dir}/data_processing/step3_post_krr_pca/output"

# Initialize input and output
if [ "$ARG1" = "normal" ] ; then
    echo "run phenotype selection with 100 rng numbers"
    phe_list="${stage_1_out_dir}/ukbb_coarse_filter_phe_list.txt"
    output="${stage_2_dir}/output"
    phe_csv="${stage_1_out_dir}/ukbb_1000_phe.csv"
else
    echo "run PCA normal selection with 100 rng numbers"
    phe_list="${stage_3_out_dir}/ukbb_pca_phe_list.txt"
    output="${stage_2_dir}/output_pca"
    phe_csv="${stage_3_out_dir}/ukbb_1000_pca_phe.csv"
fi
subj_list="${stage_1_out_dir}/ukbb_1000_subj_list.txt"
fc_mat="${stage_1_out_dir}/ukbb_1000_pfc.mat"

# Run KRR by submitting job
for j in {1..1..1}; do
    while read phe; do
        n_job="$(qstat -u ${USERNAME} | wc -l)"
        while [ $n_job -gt 50 ]; do
            echo $n_job large than 50, sleep for 1 mins
            sleep 1m
            n_job="$(qstat -u ${USERNAME} | wc -l)"
        done
        sh CBIG_MM_KRR_filter_submit_job.sh $phe $output $phe_csv $stage_2_dir $subj_list $fc_mat
    done <$phe_list
done