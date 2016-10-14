#!/bin/sh



# create MNI 1mm mask
mask_MNI1mm=GM_Mask_MNI1mm
mri_binarize --i ${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg.mgz --match 2 41 77 251 252 253 254 255 7 46 4 5 14 43 44 15 72 31 63 0 24 --inv --o $mask_MNI1mm.nii.gz


# downsample MNI 1mm mask to MNI 2mm mask
input=${mask_MNI1mm}.nii.gz
output=${mask_MNI1mm}_MNI2mm.nii.gz
MNI_ref_id=FSL_MNI152_FS4.5.0
FS_MNI_2mm=${CBIG_CODE_DIR}/data/templates/volume/FS_nonlinear_volumetric_space_4.5/gca_mean2mm.nii.gz
#FSL_MNI_2mm=${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
mri_vol2vol --mov $input --s $MNI_ref_id --targ $FS_MNI_2mm --o $output --regheader --no-save-reg




