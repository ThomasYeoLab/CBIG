#!/bin/sh
#
# There are two subjects (with two sessions of data) in the CoRR_HNU dataset involved in the example. To be able to 
# separate the data into training and test sets, these two subjects will be treated as 4 subjects (each session is 
# considered as a subject). Then the RSFC of each session will be duplicated with some random noise, resulting in 
# 8 fake subjects.
# Note that the behavioral and demographic data are faked as well.
# 
# This script will do the following things:
# 1. Generate the list of all surface data and the list of outlier files
# 2. Compute the resting-state functional connectivity matrix, based on the lists from step 1
# 3. Generate the setup file for kernel ridge regression workflow
# 4. Call the workflow wrapper to perform kernel ridge regression
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fmri_dir="$CBIG_CODE_DIR/data/example_data/CoRR_HNU"
eg_dir="$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples"
input_dir="$eg_dir/input"
##### The following two lines were used when creating the ground-truth example results
#prepare_dir="$eg_dir/output/preparation"
#KRR_dir="$eg_dir/output/KernelRidgeRegression"
################################
output_dir=$1
prepare_dir="$output_dir/preparation"
KRR_dir="$output_dir/KernelRidgeRegression"
gt_dir="$eg_dir/output/KernelRidgeRegression"
mkdir -p $prepare_dir
subjects=$(cat $input_dir/subject_list.txt)

####################################################
# Generate surf_list and motion outlier list

echo "Generating surface data list and motion outlier list ..."
lh_surf_list="$prepare_dir/lh_surf_list.txt"
rh_surf_list="$prepare_dir/rh_surf_list.txt"
outlier_list="$prepare_dir/outlier_list.txt"
if [ ! -f $lh_surf_list ] || [ ! -f $rh_surf_list ] || [ ! -f $outlier_list ]; then
	for s in $subjects; do
		orig_s=$(echo $s | cut -d'_' -f 1)
		surf_dir="$fmri_dir/$orig_s/$s/surf"
		qc_dir="$fmri_dir/$orig_s/$s/qc"
		
		run=002
		lh_surf="$surf_dir/lh.${s}_bld${run}_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"
		rh_surf="$surf_dir/rh.${s}_bld${run}_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"
		outlier="$qc_dir/${s}_bld${run}_FDRMS0.2_DVARS50_motion_outliers.txt"
		
		echo $lh_surf >> $lh_surf_list
		echo $rh_surf >> $rh_surf_list
		echo $outlier >> $outlier_list
	done
else
	echo "Warning: all the lists already exist. Skipping ..."
fi

####################################################
# Compute RSFC matrix

echo "Computing RSFC matrix ..."
ROI_dir="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/\
fsaverage6/label/"
lh_ROI="$ROI_dir/lh.Schaefer2018_400Parcels_17Networks_order.annot"
rh_ROI="$ROI_dir/rh.Schaefer2018_400Parcels_17Networks_order.annot"
RSFC_prefix="$prepare_dir/RSFC"

if [ ! -f $RSFC_prefix.mat ]; then
	# left-to-left
	cmd="matlab -nosplash -nodesktop -nodisplay -r \"CBIG_ComputeROIs2ROIsCorrelationMatrix ${RSFC_prefix}_ll.mat \
	${lh_surf_list} ${lh_surf_list} ${outlier_list} ${lh_ROI} ${lh_ROI} NONE NONE 1 0; exit;\" "
	eval $cmd
	
	# left-to-right
	cmd="matlab -nosplash -nodesktop -nodisplay -r \"CBIG_ComputeROIs2ROIsCorrelationMatrix ${RSFC_prefix}_lr.mat \
	${lh_surf_list} ${rh_surf_list} ${outlier_list} ${lh_ROI} ${rh_ROI} NONE NONE 1 0; exit;\" "
	eval $cmd
	
	# right-to-right
	cmd="matlab -nosplash -nodesktop -nodisplay -r \"CBIG_ComputeROIs2ROIsCorrelationMatrix ${RSFC_prefix}_rr.mat \
	${rh_surf_list} ${rh_surf_list} ${outlier_list} ${rh_ROI} ${rh_ROI} NONE NONE 1 0; exit;\" "
	eval $cmd

	# combine 
	# 4 subjects are not enough to do inner-loop cross-validation as well as using correlation as accuracy metric
	# Hence we fake 4 subjects by duplicating the current 4 subjects with some random noise
	cmd="matlab -nodisplay -nosplash -nodisplay -r \"
	ll = load(['${RSFC_prefix}' '_ll.mat']); \
	lr = load(['${RSFC_prefix}' '_lr.mat']); \
	rr = load(['${RSFC_prefix}' '_rr.mat']); \
	
	corr_mat = zeros(size(ll.corr_mat,1)+size(rr.corr_mat,1), size(ll.corr_mat,2)+size(rr.corr_mat,2), \
	   size(ll.corr_mat,3));
	corr_mat(1:size(ll.corr_mat, 1), 1:size(ll.corr_mat, 2), :) = ll.corr_mat; \
	corr_mat(1:size(ll.corr_mat, 1), size(ll.corr_mat, 2)+1:size(ll.corr_mat, 2)+size(rr.corr_mat, 2), :) = lr.corr_mat; \
	corr_mat(size(ll.corr_mat, 1)+1:size(ll.corr_mat, 1)+size(rr.corr_mat, 1), 1:size(ll.corr_mat, 2), :) = \
	   permute(lr.corr_mat, [2,1,3]);
	corr_mat(size(ll.corr_mat, 1)+1:size(ll.corr_mat, 1)+size(rr.corr_mat, 1), \
	   size(ll.corr_mat, 2)+1:size(ll.corr_mat, 2)+size(rr.corr_mat, 2), :) = rr.corr_mat; \
	   
	rng(1); \
	noise = rand(size(corr_mat)); \
	noise = (noise + permute(noise, [2 1 3])) ./ 2; \
	noise = bsxfun(@minus, noise, mean(noise, 3)); \
	corr_mat_fake = corr_mat + noise; \
	corr_mat = cat(3, corr_mat, corr_mat_fake); \
	
	save(['${RSFC_prefix}' '.mat'], 'corr_mat', '-v7.3'); \
	exit; \" "
	eval $cmd
	
	rm ${RSFC_prefix}_ll.mat ${RSFC_prefix}_lr.mat ${RSFC_prefix}_rr.mat
else
	echo "Warning: RSFC matrix already exists. Skipping ..."
fi

####################################################
# Prepare the setup file

fake_csv1="$input_dir/Faked_CSV1.csv"
fake_csv2="$input_dir/Faked_CSV2.csv"
fake_sub_list="$input_dir/Faked_subject_list.txt"

if [ ! -f $KRR_dir/setup_file.mat ]; then
	cmd="matlab -nosplash -nodesktop -nodisplay -r \" 
	param.sub_fold = CBIG_cross_validation_data_split( '$fake_sub_list', 'NONE', 'Subject', 'NONE', 2, 1, \
	   '$KRR_dir', ',' ); \
	
	y_names = {'Behavior_1', 'Behavior_2'}; \
	y_types = {'continuous', 'continuous'}; \
	param.y = CBIG_read_y_from_csv( {'$fake_csv1'}, 'Subject', y_names, y_types, '$fake_sub_list', \
	   fullfile('$KRR_dir', ['y.mat']), ',' ); \

	cov_names = {'Age', 'Sex'}; \
	cov_types = {'continuous', 'categorical'}; \
	param.covariates = CBIG_generate_covariates_from_csv( {'$fake_csv2'}, 'Subject', cov_names, cov_types, \
	   '$fake_sub_list', 'NONE', 'NONE', fullfile('$KRR_dir', ['covariates.mat']), ',' ); \
	
	load(['${RSFC_prefix}' '.mat']); \
	param.feature_mat = corr_mat; \
	
	param.num_inner_folds = 2; \
	param.outdir = '$KRR_dir'; \
	param.outstem = ''; \
	param.ker_param.type = 'corr'; \
	param.ker_param.scale = NaN; \
	param.lambda_set = [0 0.1 1 5 10 50 100 500 1000]; \
	param.threshold_set = []; \
	
	save(fullfile('$KRR_dir', 'setup_file.mat'), '-struct', 'param'); \
	exit; \" "
	eval $cmd
	
	rm $KRR_dir/no_relative_2_fold_sub_list.mat $KRR_dir/y.mat $KRR_dir/covariates.mat
else
	echo "Warning: setup file already exists. Skipping ..."
fi

####################################################
# Run kernel regression workflow

cmd="matlab -nodesktop -nosplash -nodisplay -r \" CBIG_KRR_workflow( fullfile('$KRR_dir', 'setup_file.mat'), 0); \
 exit; \" "
eval $cmd

####################################################
# Compare final results

echo "Comparing your results with the ground truth ..."
cmd="matlab -nosplash -nodesktop -nodisplay -r \"
your_results = load(fullfile('$KRR_dir', 'final_result.mat')); \
ground_truth = load(fullfile('$gt_dir', 'final_result.mat')); \
dif = your_results.optimal_acc - ground_truth.optimal_acc; \

if(max(abs(dif(:))) < 1e-7) \
	fprintf('Your optimal accuracies replicated the ground truth.\n');
else
	fprintf('Your optimal accuracies were different from the ground truth. Maximal difference: %f\n', max(abs(dif(:)))); \
end; \
exit; \" "
eval $cmd

