#!/bin/bash

# This script will generate example input data for 3 examples:
# 1) generate gradients 
# 2) group priors estimation
# 3) indidividual parcellation generation
# The user need to specify the output directory, which will later contain two folders:
# 1) generate gradients 
# 2) group priors estimation
# 3) indidividual parcellation generation
# In this example, we mainly focus on gMSHBM
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

##########################
# Specify output directory
##########################

out_dir=`realpath ${1}`
if [ -d $out_dir ]; then
    rm -r $out_dir
fi
mkdir -p $out_dir/estimate_group_priors
mkdir -p $out_dir/generate_individual_parcellations
mkdir -p $out_dir/generate_gradients

#########################################
# Create data lists to generate profiles
#########################################
mkdir -p $out_dir/generate_profiles_and_ini_params/data_list/fMRI_list
mkdir -p $out_dir/generate_profiles_and_ini_params/data_list/censor_list
for sess in {1..2}; do
    for sub in {1..2}; do
        # fMRI data
        lh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/surf/\
lh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $lh_fmri >> $out_dir/generate_profiles_and_ini_params/data_list/fMRI_list/lh_sub${sub}_sess${sess}.txt

        rh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/surf/\
rh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $rh_fmri >> $out_dir/generate_profiles_and_ini_params/data_list/fMRI_list/rh_sub${sub}_sess${sess}.txt

        # censor list
        censor_file="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/qc/\
subj0${sub}_sess${sess}\
_bld002_FDRMS0.2_DVARS50_motion_outliers.txt"

        echo $censor_file >> $out_dir/generate_profiles_and_ini_params/data_list/censor_list/sub${sub}_sess${sess}.txt

    done
done

#########################################
# Create data lists to generate gradients
#########################################
mkdir -p $out_dir/generate_gradients/data_list/fMRI_list
mkdir -p $out_dir/generate_gradients/data_list/censor_list
for sess in 1; do
    for sub in 1; do
        # fMRI data
        lh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/surf/\
lh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $lh_fmri >> $out_dir/generate_gradients/data_list/fMRI_list/lh_sub${sub}_sess${sess}.txt

        rh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/surf/\
rh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $rh_fmri >> $out_dir/generate_gradients/data_list/fMRI_list/rh_sub${sub}_sess${sess}.txt

        # censor list
        censor_file="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/qc/\
subj0${sub}_sess${sess}\
_bld002_FDRMS0.2_DVARS50_motion_outliers.txt"

        echo $censor_file >> $out_dir/generate_gradients/data_list/censor_list/sub${sub}_sess${sess}.txt

    done
done

##############################
# Create validation fMRI lists 
##############################
mkdir -p $out_dir/generate_individual_parcellations/data_list/validation_set/fMRI_list

for sess in 2; do
    for sub in 1; do
        # fMRI data
        lh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/\
lh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $lh_fmri >> $out_dir/generate_individual_parcellations/data_list/validation_set/fMRI_list/lh_sub${sub}.txt

        rh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/\
rh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $rh_fmri >> $out_dir/generate_individual_parcellations/data_list/validation_set/fMRI_list/rh_sub${sub}.txt

    done
done

########################################
# Generate profile lists of example data
########################################

mkdir -p $out_dir/estimate_group_priors/profile_list/training_set
mkdir -p $out_dir/generate_individual_parcellations/profile_list/test_set
mkdir -p $out_dir/generate_individual_parcellations/profile_list/validation_set

for sess in {1..2}; do
    for sub in {1..2}; do
        lh_profile="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/profiles/\
sub${sub}/sess${sess}/lh.sub${sub}_sess${sess}_fsaverage6_roifsaverage3.surf2surf_profile.nii.gz"
        echo $lh_profile >> $out_dir/estimate_group_priors/profile_list/training_set/lh_sess${sess}.txt;

        rh_profile="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/profiles/\
sub${sub}/sess${sess}/rh.sub${sub}_sess${sess}_fsaverage6_roifsaverage3.surf2surf_profile.nii.gz"
        echo $rh_profile >> $out_dir/estimate_group_priors/profile_list/training_set/rh_sess${sess}.txt;
    done
done
new_sess=0;
for sess in {1..2}; do
    for split in {1..2}; do
        ((new_sess+=1));
        for sub in {1..2}; do
            lh_profile="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/profiles/\
sub${sub}/sess${sess}/lh.sub${sub}_sess${sess}_fsaverage6_roifsaverage3.surf2surf_profile_${split}.nii.gz"
            echo $lh_profile >> $out_dir/generate_individual_parcellations/profile_list/validation_set/\
lh_sess${new_sess}.txt;
            echo $lh_profile >> $out_dir/generate_individual_parcellations/profile_list/test_set/lh_sess${new_sess}.txt;

            rh_profile="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/profiles/\
sub${sub}/sess${sess}/rh.sub${sub}_sess${sess}_fsaverage6_roifsaverage3.surf2surf_profile_${split}.nii.gz"
            echo $rh_profile >> $out_dir/generate_individual_parcellations/profile_list/validation_set/\
rh_sess${new_sess}.txt;
            echo $rh_profile >> $out_dir/generate_individual_parcellations/profile_list/test_set/rh_sess${new_sess}.txt;
        done
    done
done

########################################
# Generate gradient lists of example data
########################################

mkdir -p $out_dir/estimate_group_priors/gradient_list/training_set
mkdir -p $out_dir/generate_individual_parcellations/gradient_list/test_set
mkdir -p $out_dir/generate_individual_parcellations/gradient_list/validation_set

## we use sub1 for training set and validation set
lh_gradient="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/gradients/sub1/\
lh_emb_100_distance_matrix.mat"
echo $lh_gradient >> $out_dir/estimate_group_priors/gradient_list/training_set/gradient_list_lh.txt
echo $lh_gradient >> $out_dir/generate_individual_parcellations/gradient_list/validation_set/gradient_list_lh.txt
echo $lh_gradient >> $out_dir/generate_individual_parcellations/gradient_list/test_set/gradient_list_lh.txt

rh_gradient="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/gradients/sub1/\
rh_emb_100_distance_matrix.mat"
echo $rh_gradient >> $out_dir/estimate_group_priors/gradient_list/training_set/gradient_list_rh.txt
echo $rh_gradient >> $out_dir/generate_individual_parcellations/gradient_list/validation_set/gradient_list_rh.txt
echo $rh_gradient >> $out_dir/generate_individual_parcellations/gradient_list/test_set/gradient_list_rh.txt

## we use sub2 for training set and test set
lh_gradient="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/gradients/sub2/\
lh_emb_100_distance_matrix.mat"
echo $lh_gradient >> $out_dir/estimate_group_priors/gradient_list/training_set/gradient_list_lh.txt
echo $lh_gradient >> $out_dir/generate_individual_parcellations/gradient_list/test_set/gradient_list_lh.txt

rh_gradient="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/gradients/sub2/\
rh_emb_100_distance_matrix.mat"
echo $rh_gradient >> $out_dir/estimate_group_priors/gradient_list/training_set/gradient_list_rh.txt
echo $rh_gradient >> $out_dir/generate_individual_parcellations/gradient_list/test_set/gradient_list_rh.txt

######################################################
# Copy initialization parameters into output directory
######################################################

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/input/group" \
$out_dir/estimate_group_priors

####################################################################
# Copy estimated group priors for individual parcellation generation
####################################################################
mkdir -p $out_dir/generate_individual_parcellations/priors/gMSHBM/beta5
rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/group_priors\
/HCP_fsaverage6_40sub/400/gMSHBM/beta5/Params_Final.mat" \
$out_dir/generate_individual_parcellations/priors/gMSHBM/beta5/Params_Final.mat

####################
# Copy spatial mask
####################
mkdir -p $out_dir/generate_individual_parcellations/spatial_mask
mkdir -p $out_dir/estimate_group_priors/spatial_mask
rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask/400/\
spatial_mask_fsaverage6.mat" $out_dir/generate_individual_parcellations/spatial_mask/spatial_mask_fsaverage6.mat
rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask/400/\
spatial_mask_fsaverage6.mat" $out_dir/estimate_group_priors/spatial_mask/spatial_mask_fsaverage6.mat

echo "Example input data generated!"

