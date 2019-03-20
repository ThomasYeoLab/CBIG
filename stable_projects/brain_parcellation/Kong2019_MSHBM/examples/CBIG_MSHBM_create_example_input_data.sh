#!/bin/bash

# This script will generate example input data for two examples:
# 1) group priors estimation
# 2) indidividual parcellation generation
# The user need to specify the output directory, which will later contain two folders:
# 1) estimate_group_priors
# 2) generate_individual_parcellations
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

##########################
# Specify output directory
##########################

out_dir=$1
mkdir -p $out_dir/estimate_group_priors
mkdir -p $out_dir/generate_individual_parcellations

#########################################
# Create data lists to generate profiles
#########################################
mkdir -p $out_dir/generate_profiles_and_ini_params/data_list/fMRI_list
mkdir -p $out_dir/generate_profiles_and_ini_params/data_list/censor_list
for sess in {1,2}; do
	for sub in {1,2}; do
		# fMRI data
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/lh.subj0${sub}_sess${sess}_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz >> $out_dir/generate_profiles_and_ini_params/data_list/fMRI_list/lh_sub${sub}_sess${sess}.txt
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/rh.subj0${sub}_sess${sess}_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz >> $out_dir/generate_profiles_and_ini_params/data_list/fMRI_list/rh_sub${sub}_sess${sess}.txt
		
		# censor list
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/qc/subj0${sub}_sess${sess}_bld002_FDRMS0.2_DVARS50_motion_outliers.txt >> $out_dir/generate_profiles_and_ini_params/data_list/censor_list/sub${sub}_sess${sess}.txt
		
	done
done


######################################################
# Copy initialization parameters into output directory
######################################################

cp -r $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/examples/input/group $out_dir/estimate_group_priors

##############################
# Create validation fMRI lists 
##############################
mkdir -p $out_dir/generate_individual_parcellations/data_list/validation_fMRI_list

for sess in 2; do
	for sub in {1,2}; do
		# fMRI data
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/lh.subj0${sub}_sess${sess}_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz >> $out_dir/generate_individual_parcellations/data_list/validation_fMRI_list/lh_sub${sub}.txt
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/rh.subj0${sub}_sess${sess}_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz >> $out_dir/generate_individual_parcellations/data_list/validation_fMRI_list/rh_sub${sub}.txt
	
	done
done

########################################
# Generate profile lists of example data
########################################

mkdir -p $out_dir/estimate_group_priors/profile_list/training_set
mkdir -p $out_dir/generate_individual_parcellations/profile_list/test_set
mkdir -p $out_dir/generate_individual_parcellations/profile_list/validation_set

for sess in {1,2}; do
	for sub in {1,2}; do
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/profiles/lh.subj0${sub}.fsaverage5_roifsaverage3_BI.surf2surf_profile_cen.nii.gz >> $out_dir/estimate_group_priors/profile_list/training_set/lh_sess${sess}.txt
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/profiles/rh.subj0${sub}.fsaverage5_roifsaverage3_BI.surf2surf_profile_cen.nii.gz >> $out_dir/estimate_group_priors/profile_list/training_set/rh_sess${sess}.txt
	done
done

for split in {1,2}; do
	for sub in {1,2}; do
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess1/profiles/lh.sub${sub}_sess1_fsaverage5_roifsaverage3.surf2surf_profile_${split}.nii.gz >> $out_dir/generate_individual_parcellations/profile_list/validation_set/lh_sess${split}.txt
		echo $CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess1/profiles/rh.sub${sub}_sess1_fsaverage5_roifsaverage3.surf2surf_profile_${split}.nii.gz >> $out_dir/generate_individual_parcellations/profile_list/validation_set/rh_sess${split}.txt
	done
done

# for simplicity, the training set and test set are the same for the examples
cp -r $out_dir/estimate_group_priors/profile_list/training_set/* $out_dir/generate_individual_parcellations/profile_list/test_set

####################################################################
# Copy estimated group priors for individual parcellation generation
####################################################################

cp -r $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/examples/input/priors $out_dir/generate_individual_parcellations

	
echo "Example input data generated!"

