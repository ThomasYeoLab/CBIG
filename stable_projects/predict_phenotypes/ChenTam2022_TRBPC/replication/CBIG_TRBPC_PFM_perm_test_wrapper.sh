#!/bin/bash
# this function runs all the kernel regressions in Chen & Tam 2021 paper
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

pred_results_dir=$1
perm_out_dir=$2
###################################################################################################
# set common variables
###################################################################################################

scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
repdata_dir=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
y_list=${repdata_dir}/KRR/variables_to_predict.txt
site_list=${repdata_dir}/KRR/sites_subjects_all_task_all_score_unrelated.txt
feature_file=${repdata_dir}/KRR/feature_files_allFC.mat
outstem=all_score
cluster=CBIG_cluster
num_behav=`cat $y_list | wc -l`
perm_num=200
total_num=100000

###################################################################################################
# run permutation test for PFM of multi-kernel regression
###################################################################################################
$scripts_dir/permutation/PFM/CBIG_TRBPC_PFM_perm_wrapper.sh -i ${pred_results_dir} -f ${feature_file} -s ${outstem} \
-n ${perm_num} -t ${total_num} -o ${perm_out_dir} -g ${site_list} -c ${cluster}
