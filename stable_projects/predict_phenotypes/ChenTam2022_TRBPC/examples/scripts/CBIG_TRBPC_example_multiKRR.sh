#!/bin/bash
# this function runs all the example multi-kernel regression
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# set cluster; change it to cluster="none" if you don't have a cluster
cluster=CBIG_cluster

# set other input variables
outdir=$1/multiKRR
mkdir -p $outdir
LF=$outdir/multiKRR.log

scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
input_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/examples/input
subject_list=${input_dir}/subject_list.txt
y_list=${input_dir}/y_list.txt
FD_file=none
DVARS_file=none
outstem=all_score
csv_file=${input_dir}/behav_covariates.csv
feature_file=${outdir}/FC_all.mat
covariate_list=${input_dir}/covariates.txt

# run single kernel KRR
matlab -nodesktop -nosplash -nodisplay -r " addpath(genpath('$scripts_dir')); CBIG_TRBPC_example_multiKRR_get_inputs(\
    '${outdir}');exit; " >> $LF 2>&1

${scripts_dir}/KRR_LpOCV/MultiKernel/CBIG_TRBPC_multiKRR_LpOCV_workflow.sh -subject_list ${subject_list} \
-feature_files ${feature_file} -outstem ${outstem} -y_list ${y_list} -covariate_list ${covariate_list} \
-FD_file ${FD_file} -DVARS_file ${DVARS_file} -outdir ${outdir} -csv_file ${csv_file} -cluster "$cluster"
