#!/bin/bash

# This script will generate replication input data for:
# 1) group priors estimation
# 2) individual parcellation generation
# The user need to specify the output directory, which will later contain two folders:
# 1) estimate_group_priors
# 2) generate_individual_parcellations
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

##########################
# Specify output directory
##########################

out_dir=`realpath ${1}`
echo $out_dir
mkdir -p $out_dir/estimate_group_priors/GSP_HNU
mkdir -p $out_dir/estimate_group_priors/HCP
mkdir -p $out_dir/generate_individual_parcellations/GSP_HNU
mkdir -p $out_dir/generate_individual_parcellations/HCP


######################################################
# Copy initialization parameters into output directory
######################################################

cp -r "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input/\
estimate_group_priors/GSP_HNU/group" $out_dir/estimate_group_priors/GSP_HNU

cp -r "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input/\
estimate_group_priors/HCP/group" $out_dir/estimate_group_priors/HCP

########################################
# Generate profile lists of training data
########################################

mkdir -p $out_dir/estimate_group_priors/GSP_HNU/profile_list/training_set
mkdir -p $out_dir/estimate_group_priors/HCP/profile_list/training_set
mkdir -p $out_dir/generate_individual_parcellations/GSP_HNU/profile_list/test_set
mkdir -p $out_dir/generate_individual_parcellations/HCP/profile_list/test_set

### Generate profile lists of GSP training set
GSP_subject_list="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input/\
estimate_group_priors/GSP_HNU/subject_list/GSP_37subname_training_set_list.txt"

for sess in {1,2}; do
	lh_GSP_train_profile_list=$out_dir/estimate_group_priors/GSP_HNU/profile_list/training_set/lh_sess${sess}.txt
	rh_GSP_train_profile_list=$out_dir/estimate_group_priors/GSP_HNU/profile_list/training_set/rh_sess${sess}.txt
	
	cat $GSP_subject_list | while read sub; do	
		cd $CBIG_MSHBM_REP_GSP_DIR/$sub/${sub}_Ses$sess/ruby_surf2surf_profiles; 
		for i in lh.$sub.fsaverage5_roifsaverage3_BI.surf2surf_profile.nii.gz; do
			echo $CBIG_MSHBM_REP_GSP_DIR/$sub/${sub}_Ses$sess/ruby_surf2surf_profiles/$i >> $lh_GSP_train_profile_list;
		done;
		
		for j in rh.$sub.fsaverage5_roifsaverage3_BI.surf2surf_profile.nii.gz; do
			echo $CBIG_MSHBM_REP_GSP_DIR/$sub/${sub}_Ses$sess/ruby_surf2surf_profiles/$j >> $rh_GSP_train_profile_list;
		done;
	done;
done;

### Generate profile lists of HCP training set
train_subject_list="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input/\
estimate_group_priors/HCP/subject_list/HCP_40subname_training_set_list.txt"

cat $train_subject_list | while read sub; do
	HCP_train_profile_list=$out_dir/estimate_group_priors/HCP/profile_list/training_set;
	
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_LR_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess1.txt;
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST1_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_RL_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess2.txt;
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST2_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST2_LR_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess3.txt;
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST2_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST2_RL_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess4.txt;
	
done


####################################################################
# Copy estimated group priors for individual parcellation generation
####################################################################

cp -r "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input\
/generate_individual_parcellations/GSP_HNU/priors" $out_dir/generate_individual_parcellations/GSP_HNU

cp -r "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input\
/generate_individual_parcellations/HCP/priors" $out_dir/generate_individual_parcellations/HCP

########################################
# Generate profile lists of test data
########################################

### Generate profile lists of HNU test set

for sess in {1..10}; do
	lh_HNU_train_profile_list="$out_dir/generate_individual_parcellations/GSP_HNU/profile_list/test_set/\
lh_sess${sess}.txt"
	rh_HNU_train_profile_list="$out_dir/generate_individual_parcellations/GSP_HNU/profile_list/test_set/\
rh_sess${sess}.txt"
	
	for sub in 0{1..9} {10..30}; do
	
		cd $CBIG_MSHBM_REP_HNU_DIR/subj$sub/subj${sub}_sess$sess/ruby_surf2surf_profiles; 
		for i in lh.subj${sub}.fsaverage5_roifsaverage3_BI.surf2surf_profile_cen.nii.gz; do
			echo "$CBIG_MSHBM_REP_HNU_DIR/subj$sub/subj${sub}_sess$sess/ruby_surf2surf_profiles/$i" \
>> $lh_HNU_train_profile_list;
		done;

		for j in rh.subj${sub}.fsaverage5_roifsaverage3_BI.surf2surf_profile_cen.nii.gz; do
			echo "$CBIG_MSHBM_REP_HNU_DIR/subj$sub/subj${sub}_sess$sess/ruby_surf2surf_profiles/$j" \
>> $rh_HNU_train_profile_list;
		done;
	done;
done;

### Generate profile lists of HCP test set
test_subject_list="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication/input/\
generate_individual_parcellations/HCP/subject_list/HCP_755subname_test_set_list.txt"

cat $test_subject_list | while read sub; do
	HCP_test_profile_list=$out_dir/generate_individual_parcellations/HCP/profile_list/test_set;
	
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_LR_downsampleROI_profile.mat" \
>> $HCP_test_profile_list/sess1.txt;
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST1_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_RL_downsampleROI_profile.mat" \
>> $HCP_test_profile_list/sess2.txt;
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST2_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST2_LR_downsampleROI_profile.mat" \
>> $HCP_test_profile_list/sess3.txt;
	echo "$CBIG_MSHBM_REP_HCP_DIR/$sub/MNINonLinear/Results/\
rfMRI_REST2_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST2_RL_downsampleROI_profile.mat" \
>> $HCP_test_profile_list/sess4.txt;
	
done

	
echo "Replication input data generated!"

