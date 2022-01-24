#!/bin/bash
# this function computes the PFM
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# set cluster; change it to cluster="none" if you don't have a cluster
cluster=CBIG_cluster

###################################################################################################
# set common variables
###################################################################################################
results_dir=$1
scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
input_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/examples/input
y_list=${input_dir}/y_list.txt
cluster=circ-spool
num_behav=`cat $y_list | wc -l`

###################################################################################################
# compute PFM for multi-kernel regression
###################################################################################################
model=multiKRR
krr_dir=$results_dir/multiKRR
feature_file=${krr_dir}/FC_all.mat
sub_fold=$results_dir/multiKRR/no_relative_3_fold_sub_list.mat
outdir=$results_dir/PFM/multiKRR
$scripts_dir/PFM/CBIG_TRBPC_PFM_wrapper.sh -i ${krr_dir} -f ${feature_file} -s ${sub_fold} -n ${num_behav} \
-o ${outdir} -m ${model} -c ${cluster}

###################################################################################################
# compute PFM for single-kernel regression
###################################################################################################
model=singleKRR
krr_dir=$results_dir/singleKRR
feature_file=${input_dir}/FC_rs.mat
sub_fold=$results_dir/singleKRR/no_relative_3_fold_sub_list.mat
outdir=$results_dir/PFM/singleKRR
$scripts_dir/PFM/CBIG_TRBPC_PFM_wrapper.sh -i ${krr_dir} -f ${feature_file} -s ${sub_fold} -n ${num_behav} \
-o ${outdir} -m ${model} -c ${cluster}

###################################################################################################
# compute PFM for single-kernel regression
###################################################################################################

model=LRR
lrr_dir=$results_dir/LRR
feature_file=${input_dir}/FC_rs.mat
sub_fold=$results_dir/LRR/no_relative_3_fold_sub_list.mat
outdir=$results_dir/PFM/LRR
$scripts_dir/PFM/CBIG_TRBPC_PFM_wrapper.sh -i ${lrr_dir} -f ${feature_file} -s ${sub_fold} -n ${num_behav} \
-o ${outdir} -m ${model} -c ${cluster}
