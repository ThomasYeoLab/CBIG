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
mkdir -p $out_dir/estimate_group_priors/HCP_fs_LR_32k
mkdir -p $out_dir/estimate_group_priors/HCP_fsaverage6
mkdir -p $out_dir/generate_individual_parcellations/HCP_fs_LR_32k
mkdir -p $out_dir/generate_individual_parcellations/ABCD_fsaverage6


######################################################
# Copy initialization parameters into output directory
######################################################

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/input/\
estimate_group_priors/HCP_fs_LR_32k/group" $out_dir/estimate_group_priors/HCP_fs_LR_32k

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/input/\
estimate_group_priors/HCP_fsaverage6/group" $out_dir/estimate_group_priors/HCP_fsaverage6

####################
# Copy spatial mask
####################

mkdir -p $out_dir/estimate_group_priors/HCP_fs_LR_32k/spatial_mask
mkdir -p $out_dir/estimate_group_priors/HCP_fsaverage6/spatial_mask

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask/400/\
spatial_mask_fs_LR_32k.mat" $out_dir/estimate_group_priors/HCP_fs_LR_32k/spatial_mask/spatial_mask_fs_LR_32k.mat

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask/400/\
spatial_mask_fsaverage6.mat" $out_dir/estimate_group_priors/HCP_fsaverage6/spatial_mask/spatial_mask_fsaverage6.mat


######################################################
# Generate profile and gradient lists of training data
######################################################

mkdir -p $out_dir/estimate_group_priors/HCP_fs_LR_32k/profile_list/training_set
mkdir -p $out_dir/estimate_group_priors/HCP_fsaverage6/profile_list/training_set
mkdir -p $out_dir/generate_individual_parcellations/HCP_fs_LR_32k/profile_list/test_set
mkdir -p $out_dir/generate_individual_parcellations/ABCD_fsaverage6/profile_list/test_set
mkdir -p $out_dir/estimate_group_priors/HCP_fs_LR_32k/gradient_list/training_set
mkdir -p $out_dir/estimate_group_priors/HCP_fsaverage6/gradient_list/training_set
mkdir -p $out_dir/generate_individual_parcellations/HCP_fs_LR_32k/gradient_list/test_set
mkdir -p $out_dir/generate_individual_parcellations/ABCD_fsaverage6/gradient_list/test_set

### Generate profile lists of HCP training set
train_subject_list="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/input/\
estimate_group_priors/subject_list/HCP_40subname_training_set_list.txt"

cat $train_subject_list | while read sub; do
    HCP_train_profile_list=$out_dir/estimate_group_priors/HCP_fs_LR_32k/profile_list/training_set;

    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_LR_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess1.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST1_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_RL_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess2.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST2_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST2_LR_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess3.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST2_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST2_RL_downsampleROI_profile.mat" \
>> $HCP_train_profile_list/sess4.txt;

done

cat $train_subject_list | while read sub; do
    HCP_train_profile_list=$out_dir/estimate_group_priors/HCP_fsaverage6/profile_list/training_set;

    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/\
lh.rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/lh_sess1.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST1_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/\
lh.rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/lh_sess2.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST2_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/\
lh.rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/lh_sess3.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST2_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/\
lh.rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/lh_sess4.txt;

    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/\
rh.rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/rh_sess1.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST1_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/\
rh.rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/rh_sess2.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST2_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/\
rh.rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/rh_sess3.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/S1200/individuals/$sub/MNINonLinear/Results/\
rfMRI_REST2_RL/postprocessing/MSM_reg_wbsgrayordinatecortex/\
rh.rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean_regress_fs6_fs3_profile.mat" \
>> $HCP_train_profile_list/rh_sess4.txt;

done

### Generqate gradient lists of HCP training set
cat $train_subject_list | while read sub; do
    HCP_train_gradient_list=$out_dir/estimate_group_priors/HCP_fs_LR_32k/gradient_list/training_set;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/GradientMaps_Gordon/Speed_Up_Version/HCP_dataA_separate_runs/\
0_10_200_3.2_3_10/$sub/lh_emb_1000_distance_matrix.mat" >> $HCP_train_gradient_list/gradient_list_lh.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/GradientMaps_Gordon/Speed_Up_Version/HCP_dataA_separate_runs/\
0_10_200_3.2_3_10/$sub/rh_emb_1000_distance_matrix.mat" >> $HCP_train_gradient_list/gradient_list_rh.txt;

done

cat $train_subject_list | while read sub; do
    HCP_train_gradient_list=$out_dir/estimate_group_priors/HCP_fsaverage6/gradient_list/training_set;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/GradientMaps_Gordon/CBIG_Speed_Up_Version_FS6/100_200_3.2/\
HCP_dataA/$sub/lh_emb_100_distance_matrix.mat" >> $HCP_train_gradient_list/gradient_list_lh.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/GradientMaps_Gordon/CBIG_Speed_Up_Version_FS6/100_200_3.2/\
HCP_dataA/$sub/rh_emb_100_distance_matrix.mat" >> $HCP_train_gradient_list/gradient_list_rh.txt;

done


####################
# Copy spatial mask
####################

mkdir -p $out_dir/generate_individual_parcellations/HCP_fs_LR_32k/spatial_mask
mkdir -p $out_dir/generate_individual_parcellations/ABCD_fsaverage6/spatial_mask

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask/400/\
spatial_mask_fs_LR_32k.mat" $out_dir/generate_individual_parcellations/HCP_fs_LR_32k/spatial_mask/\
spatial_mask_fs_LR_32k.mat

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask/400/\
spatial_mask_fsaverage6.mat" $out_dir/generate_individual_parcellations/ABCD_fsaverage6/spatial_mask/\
spatial_mask_fsaverage6.mat


####################################################################
# Copy estimated group priors for individual parcellation generation
####################################################################
mkdir -p $out_dir/generate_individual_parcellations/HCP_fs_LR_32k/priors/gMSHBM/beta50
mkdir -p $out_dir/generate_individual_parcellations/ABCD_fsaverage6/priors/gMSHBM/beta5

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/group_priors/HCP_fs_LR_32k_40sub/\
400/gMSHBM/beta50/Params_Final.mat" $out_dir/generate_individual_parcellations/HCP_fs_LR_32k/priors/gMSHBM/beta50/\
Params_Final.mat

rsync -az "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/group_priors/HCP_fsaverage6_40sub/\
400/gMSHBM/beta5/Params_Final.mat" $out_dir/generate_individual_parcellations/ABCD_fsaverage6/priors/gMSHBM/beta5/\
Params_Final.mat

##################################################
# Generate profile and gradient lists of test data
##################################################

### Generate profile lists of HCP test set
test_subject_list="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/input/\
estimate_group_priors/subject_list/HCP_755subname_test_set_list.txt"
for sess in {1..12}; do
    cat $test_subject_list | while read sub; do
        HCP_test_profile_list=$out_dir/generate_individual_parcellations/HCP_fs_LR_32k/profile_list/test_set;
        echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/HCP_fake_multisession_profiles/\
3/${sub}/${sub}_profile_multisess${sess}.mat" >> $HCP_test_profile_list/sess${sess}.txt;
    done
done

### Generqate gradient lists of HCP test set
cat $test_subject_list | while read sub; do
    HCP_test_gradient_list=$out_dir/generate_individual_parcellations/HCP_fs_LR_32k/gradient_list/test_set;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/GradientMaps_Gordon/Speed_Up_Version/HCP_dataC_separate_runs/\
0_10_200_3.2_3_10/$sub/lh_emb_1000_distance_matrix.mat" >> $HCP_test_gradient_list/gradient_list_lh.txt;
    echo "$CBIG_ArealMSHBM_REP_HCP_DIR/HCP_derivatives/GradientMaps_Gordon/Speed_Up_Version/HCP_dataC_separate_runs/\
0_10_200_3.2_3_10/$sub/rh_emb_1000_distance_matrix.mat" >> $HCP_test_gradient_list/gradient_list_rh.txt;

done

### Generate profile lists of ABCD test set
test_subject_list="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/input/\
estimate_group_priors/subject_list/ABCD_5subname_test_set_list.txt"

for sess in {1..4}; do
    cat $test_subject_list | while read sub; do
        ABCD_test_profile_list=$out_dir/generate_individual_parcellations/ABCD_fsaverage6/profile_list/test_set;
        
        echo "$CBIG_REPDATA_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/ABCD_profile/$sub\
/lh.bld00${sess}_sm6_fs6_fs3_profile.mat" >> $ABCD_test_profile_list/lh_sess${sess}.txt
        echo "$CBIG_REPDATA_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/ABCD_profile/$sub\
/rh.bld00${sess}_sm6_fs6_fs3_profile.mat" >> $ABCD_test_profile_list/rh_sess${sess}.txt
    done
done

### Generqate gradient lists of ABCD test set
cat $test_subject_list | while read sub; do
    ABCD_test_gradient_list=$out_dir/generate_individual_parcellations/ABCD_fsaverage6/gradient_list/test_set;
    echo "$CBIG_REPDATA_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/CBIG_Speed_Up_Version_FS6/\
100_200_3.2/ABCD/$sub/lh_emb_100_distance_matrix.mat" >> $ABCD_test_gradient_list/gradient_list_lh.txt;
    echo "$CBIG_REPDATA_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/CBIG_Speed_Up_Version_FS6/\
100_200_3.2/ABCD/$sub/rh_emb_100_distance_matrix.mat" >> $ABCD_test_gradient_list/gradient_list_rh.txt;

done


echo "Replication input data generated!"

