#!/bin/bash
# this function runs all the control analysis of regression models
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

top_outdir=$1
###################################################################################################
# set common variables
###################################################################################################

scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
repdata_dir=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
subject_list=${repdata_dir}/lists/release2_subjects_all_task_all_score_unrelated.txt
y_list=${repdata_dir}/KRR/variables_to_predict.txt
FD_file=none
DVARS_file=none
outstem=all_score
csv_file=${repdata_dir}/KRR/score_motion_all_subjects.csv

###################################################################################################
# run single-kernel regression with age and sex as covariates
###################################################################################################

for curr_task in rs mid nback sst meanFC
do
	feature_file=${repdata_dir}/FC/FC_subjects_all_task_all_score_${curr_task}_bp_sm6.mat
	covariate_list=${repdata_dir}/KRR/covariates_${curr_task}_age_sex.txt
	outdir=${top_outdir}/KRR_regress_age_sex/${curr_task}
	${scripts_dir}/KRR_LpOCV/SingleKernel/CBIG_TRBPC_KRR_LpOCV_workflow.sh -subject_list ${subject_list} \
	-feature_file ${feature_file} -outstem ${outstem} -y_list ${y_list} -covariate_list ${covariate_list} \
	-FD_file ${FD_file} -DVARS_file ${DVARS_file} -outdir ${outdir} -csv_file ${csv_file}
done

###################################################################################################
# run multi-kernel regression with age and sex as covariates
###################################################################################################
feature_files=$repdata_dir/KRR/feature_files_allFC.mat
covariate_list=${repdata_dir}/KRR/covariates_all_age_sex.txt
outdir=${top_outdir}/KRR_regress_age_sex/allFC

${scripts_dir}/KRR_LpOCV/MultiKernel/CBIG_TRBPC_multiKRR_LpOCV_workflow.sh -subject_list ${subject_list} \
-feature_files ${feature_files} -outstem ${outstem} -y_list ${y_list} -covariate_list ${covariate_list} -FD_file \
${FD_file} -DVARS_file ${DVARS_file} -outdir ${outdir} -csv_file ${csv_file}

###################################################################################################
# run linear ridge regression
###################################################################################################
for curr_task in rs mid nback sst allFC
do
	feature_file=${repdata_dir}/FC/FC_norm_${curr_task}.mat
	covariate_list=${repdata_dir}/KRR/covariates_${curr_task}.txt
	outdir=${top_outdir}/LRR/${curr_task}
	domain=${repdata_dir}/LRR/domain.mat
	${scripts_dir}/LRR_LpOCV/CBIG_TRBPC_LRR_LpOCV_workflow.sh -subject_list ${subject_list} \
	-feature_file ${feature_file} -outstem ${outstem} -y_list ${y_list} -covariate_list ${covariate_list} \
	-FD_file ${FD_file} -DVARS_file ${DVARS_file} -outdir ${outdir} -csv_file ${csv_file} -domain ${domain}
done