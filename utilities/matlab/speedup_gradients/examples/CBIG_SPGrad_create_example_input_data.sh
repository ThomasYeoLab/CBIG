#!/bin/bash

# This script will generate example input data
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

##########################
# Specify output directory
##########################

out_dir=`realpath ${1}`
if [ -d $out_dir ]; then
    rm -r $out_dir
fi


#########################################
# Create data lists to generate gradientss
#########################################
mkdir -p $out_dir/data_list/fMRI_list
mkdir -p $out_dir/data_list/censor_list
for sess in {1,2}; do
    for sub in 1; do
        # fMRI data
        lh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/surf/\
lh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $lh_fmri >> $out_dir/data_list/fMRI_list/lh_sub${sub}.txt

        rh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/surf/\
rh.subj0${sub}_sess${sess}\
_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz"

        echo $rh_fmri >> $out_dir/data_list/fMRI_list/rh_sub${sub}.txt

        # censor list
        censor_file="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/\
subj0$sub/subj0${sub}_sess${sess}/qc/\
subj0${sub}_sess${sess}\
_bld002_FDRMS0.2_DVARS50_motion_outliers.txt"

        echo $censor_file >> $out_dir/data_list/censor_list/sub${sub}.txt

    done
done

echo "Example input data generated!"

