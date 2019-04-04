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
# 3. Call $CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/utilities/
#    CBIG_LiGSR_compute_FSM_from_FC.m to compute functional similarity matrix across subjects.
# 4. Call $CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/utilities/
#    CBIG_LiGSR_NchooseD_families.m to randomly generate 10 delete-2 jackknife samples.
# 5. Call 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSRexamples/scripts/CBIG_LiGSR_example_explained_variance.m
#    to estimate the explained variance on the full set as well as on each jackknife sample.
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fmri_dir="$CBIG_CODE_DIR/data/example_data/CoRR_HNU"
eg_dir="$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples"
input_dir="$eg_dir/input"
##### The following two lines were used when creating the ground-truth example results
#prepare_dir="$eg_dir/output/preparation"
#LME_dir="$eg_dir/output/VarianceComponentModel"
################################
output_dir=$1
prepare_dir="$output_dir/preparation"
LME_dir="$output_dir/VarianceComponentModel"
gt_dir="$eg_dir/output/VarianceComponentModel"
mkdir -p $prepare_dir
subjects=$(cat $input_dir/subject_list.txt)

fake_csv1="$input_dir/Faked_CSV1.csv"
fake_csv2="$input_dir/Faked_CSV2.csv"
fake_sub_list="$input_dir/Faked_subject_list.txt"

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
	   size(ll.corr_mat,3));\
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

########################################
# Compute FSM from RSFC
echo "Computing functional similarity matrix (FSM) ..."
FSM_file="$LME_dir/FSM.mat";
if [ ! -f $FSM_file ]; then
	cmd="matlab -nodesktop -nosplash -nodisplay -r \" addpath(fullfile('$CBIG_CODE_DIR', 'stable_projects', \
	   'preprocessing', 'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities')); \
	   CBIG_LiGSR_compute_FSM_from_FC('$RSFC_prefix.mat', '$FSM_file', 'corr'); exit; \" "
	eval $cmd
else
	echo "Warning: FSM already exists. Skipping ..."
fi


########################################
# Estimate explained trait variance on the full set
echo "Estimating explained variance on the full set ..."
cmd="matlab -nodesktop -nosplash -nodisplay -r \" addpath(fullfile('$CBIG_CODE_DIR', 'stable_projects', \
   'preprocessing', 'Li2019_GSR', 'examples', 'scripts')); CBIG_LiGSR_example_explained_variance( '$fake_csv1', \
   '$fake_csv2', '$fake_sub_list', 'NONE', '$FSM_file', '$LME_dir', 'fullset' ) ; exit;\" "
eval $cmd

########################################
# Generate jackknife samples (delete 2 subjects, repeated 10 times)
echo "Generating delete-d jackknife sample ..."
cmd="matlab -nodesktop -nosplash -nodisplay -r \" addpath(fullfile('$CBIG_CODE_DIR', 'stable_projects', \
   'preprocessing', 'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities')); CBIG_LiGSR_NchooseD_families( \
   '$fake_sub_list', 'none', 'none', 'none', 2, 10, '$LME_dir/jackknife_lists', 'Toy8subjects' ); exit; \" "
eval $cmd

########################################
# Estimate explained variance for each jakknife sample
echo "Estimating explained variance on each jackknife sample ..."
cmd="matlab -nodesktop -nosplash -nodisplay -r \" 
   addpath(fullfile('$CBIG_CODE_DIR', 'stable_projects', 'preprocessing', 'Li2019_GSR', 'examples', 'scripts'));
   for i = 1:10
       rm_sub_list = fullfile('$LME_dir', 'jackknife_lists', ['Toy8subjects_choose2_set' num2str(i) '.txt']); \
       CBIG_LiGSR_example_explained_variance( '$fake_csv1', '$fake_csv2', '$fake_sub_list', \
          rm_sub_list, '$FSM_file', '$LME_dir', ['del2_set' num2str(i)] );
   end; \
   exit; \" "
eval $cmd

#########################################
# Compare your results with ground truth
echo "Comparing your results with the ground truth ..."
cmd="matlab -nodesktop -nosplash -nodisplay -r \" 
for n = 1:2 
	yours = load(fullfile('$LME_dir', 'fullset', sprintf('m2_QuantileNorm_Behavior_%d.mat', n))); 
	gt = load(fullfile('$gt_dir', 'fullset', sprintf('m2_QuantileNorm_Behavior_%d.mat', n))); 
	
	if((yours.morpho.m2 - gt.morpho.m2) < 1e-15) || (isnan(yours.morpho.m2) && isnan(gt.morpho.m2))
		fprintf('Your estimation of explained variance of behavior %d using the full set was correct.\n', n); 
	else 
		fprintf('Your estimation of explained variance of behavior %d using the full set was different from the \
ground truth by %f.\n', n, yours.morpho.m2 - gt.morpho.m2); 
	end; 
	
	if((yours.morpho.SE - gt.morpho.SE) < 1e-15) || (isnan(yours.morpho.SE) && isnan(gt.morpho.SE))
		fprintf('Your estimation of standard deviation for behavior %d using the full set was correct.\n', n); \
	else 
		fprintf('Your estimation of standard deviation for behavior %d using the full set was different from the \
ground truth by %f.\n', n, yours.morpho.SE - gt.morpho.SE); 
	end; 
end; 

for set = 1:10
	for n = 1:2 
		yours = load(fullfile('$LME_dir', sprintf('del2_set%d', set), sprintf('m2_QuantileNorm_Behavior_%d.mat', n))); 
		gt = load(fullfile('$gt_dir', sprintf('del2_set%d', set), sprintf('m2_QuantileNorm_Behavior_%d.mat', n))); 
		
		if((yours.morpho.m2 - gt.morpho.m2) < 1e-15) || (isnan(yours.morpho.m2) && isnan(gt.morpho.m2))
			fprintf('Your estimation of explained variance of behavior %d using the jackknife sample %d was \
correct.\n', n, set); 
		else 
			fprintf('Your estimation of explained variance of behavior %d using the jackknife sample %d was \
different from the ground truth by %f.\n', n, set, yours.morpho.m2 - gt.morpho.m2); \
		end; 
		
		if((yours.morpho.SE - gt.morpho.SE) < 1e-15) || (isnan(yours.morpho.SE) && isnan(gt.morpho.SE))
			fprintf('Your estimation of standard deviation for behavior %d using the jackknife sample %d was \
correct.\n', n, set); 
		else 
			fprintf('Your estimation of standard deviation for behavior %d using the jackknife sample %d was \
different from the ground truth by %f.\n', n, set, yours.morpho.SE - gt.morpho.SE); 
		end; 
	end; 
end; \

exit; \" "
eval $cmd


