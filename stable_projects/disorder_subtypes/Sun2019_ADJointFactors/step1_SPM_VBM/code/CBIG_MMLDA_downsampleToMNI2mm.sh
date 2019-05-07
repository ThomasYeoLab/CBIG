#!/bin/sh

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

input=$1

ext=${input#*.}
filename=${input%%.*}
output=${filename}_MNI2mm.${ext}
MNI_ref_id=FSL_MNI152_FS4.5.0
MNI_temp_2mm=${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz

mri_vol2vol --mov $input --s $MNI_ref_id --targ $MNI_temp_2mm --o $output --regheader --no-save-reg
