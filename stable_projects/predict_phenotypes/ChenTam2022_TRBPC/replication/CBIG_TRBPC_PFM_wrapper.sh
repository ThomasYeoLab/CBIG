#!/bin/bash
# this function computes the PFM for all models used in our paper
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

pred_results_dir=$1
PFM_outdir=$2
###################################################################################################
# set common variables
###################################################################################################

scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
repdata_dir=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
y_list=${repdata_dir}/KRR/variables_to_predict.txt
cluster=CBIG_cluster
num_behav=`cat $y_list | wc -l`
N_per_job=1000
total_num=2000

###################################################################################################
# compute PFM for multi-kernel regression
###################################################################################################
model=multiKRR
for curr_task in allFC
do
	input_dir=$pred_results_dir/KRR/$curr_task
    	feature_file=${repdata_dir}/KRR/feature_files_$curr_task.mat
    	sub_fold=$pred_results_dir/KRR/$curr_task/no_relative_3_fold_sub_list.mat
	outdir=$PFM_outdir/KRR/$curr_task
	$scripts_dir/PFM/CBIG_TRBPC_PFM_wrapper.sh -i ${input_dir} -f ${feature_file} -s ${sub_fold} -n ${num_behav} \
    -o ${outdir} -m ${model} -c ${cluster}
done

###################################################################################################
# compute PFM for single-kernel regression
###################################################################################################
model=singleKRR
for curr_task in rs mid nback sst
do
	input_dir=$pred_results_dir/KRR/$curr_task
    	feature_file=${repdata_dir}/FC/FC_subjects_all_task_all_score_${curr_task}_bp_sm6.mat
    	sub_fold=$pred_results_dir/KRR/$curr_task/no_relative_3_fold_sub_list.mat
	outdir=$PFM_outdir/KRR/$curr_task
	$scripts_dir/PFM/CBIG_TRBPC_PFM_wrapper.sh -i ${input_dir} -f ${feature_file} -s ${sub_fold} -n ${num_behav} \
    -o ${outdir} -m ${model} -c ${cluster}
done

###################################################################################################
# compute PFM for single-kernel regression
###################################################################################################
model=LRR
for curr_task in rs mid nback sst allFC
do
	input_dir=$pred_results_dir/LRR/$curr_task
    	feature_file=${repdata_dir}/FC/FC_norm_${curr_task}.mat
    	sub_fold=$pred_results_dir/LRR/$curr_task/no_relative_3_fold_sub_list.mat
	outdir=$PFM_outdir/LRR/$curr_task
	$scripts_dir/PFM/CBIG_TRBPC_PFM_wrapper.sh -i ${input_dir} -f ${feature_file} -s ${sub_fold} -n ${num_behav} \
    -o ${outdir} -m ${model} -c ${cluster}
done
