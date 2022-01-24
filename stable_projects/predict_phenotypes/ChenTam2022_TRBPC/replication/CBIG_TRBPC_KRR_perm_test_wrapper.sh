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
outstem=all_score
cluster=CBIG_cluster
num_behav=`cat $y_list | wc -l`
N_per_job=1000
total_num=2000
###################################################################################################
# run permutation test for single-kernel regression
###################################################################################################
model=singleKRR
for curr_task in rs mid nback sst
do
	input_dir=$pred_results_dir/$curr_task
	outdir=$perm_out_dir/$curr_task
	$scripts_dir/permutation/KRR/CBIG_TRBPC_KRR_perm_wrapper.sh -i ${input_dir} -b ${num_behav} -s ${outstem} \
	-n ${N_per_job} -t ${total_num} -o ${outdir} -m ${model} -g ${site_list} -c ${cluster}
done

###################################################################################################
# run multi-kernel regression
###################################################################################################
model=multiKRR
for curr_task in allFC
do
	input_dir=$pred_results_dir/$curr_task
	outdir=$perm_out_dir/$curr_task
	$scripts_dir/permutation/KRR/CBIG_TRBPC_KRR_perm_wrapper.sh -i ${input_dir} -b ${num_behav} -s ${outstem} \
	-n ${N_per_job} -t ${total_num} -o ${outdir} -m ${model} -g ${site_list} -c ${cluster}
done
