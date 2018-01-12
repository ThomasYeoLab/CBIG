#!/bin/sh
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

curr_SUBJECTS_DIR=$SUBJECTS_DIR
export SUBJECTS_DIR=$CBIG_CODE_DIR/data/templates/volume

cmd="mri_vol2vol --mov $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg.mgz --s FSL_MNI152_FS4.5.0 --targ $CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_1mm_brain.nii.gz --o $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg_182x218x182.nii.gz --regheader --no-save-reg"
echo $cmd
eval $cmd

export SUBJECTS_DIR=$curr_SUBJECTS_DIR
