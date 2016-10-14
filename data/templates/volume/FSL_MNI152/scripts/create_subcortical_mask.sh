#!/bin/sh



# create MNI 1mm mask
mask_MNI1mm=SubcorticalLooseMask_MNI1mm
mri_binarize --i ${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg.mgz --match 16 --match 8 --match 9 --match 10 --match 11 --match 12 --match 13 --match 17 --match 18 --match 26 --match 27 --match 28 --match 47 --match 48 --match 49 --match 50 --match 51 --match 52 --match 53 --match 54 --match 58 --match 59 --match 60 --o $mask_MNI1mm.nii.gz

# smooth MNI 1mm mask
MNI_sm=6
mask_MNI1mm_sm=${mask_MNI1mm}_sm${MNI_sm}
std=`echo ${MNI_sm}/2/2.35482 | bc -l`
echo $std
fslmaths $mask_MNI1mm.nii.gz -s $std $mask_MNI1mm_sm.nii.gz

# downsample MNI 1mm mask to MNI 2mm mask
input=${mask_MNI1mm_sm}.nii.gz
output=${mask_MNI1mm_sm}_MNI2mm.nii.gz
MNI_ref_id=FSL_MNI152_FS4.5.0
MNI_temp_2mm=${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
mri_vol2vol --mov $input --s $MNI_ref_id --targ $MNI_temp_2mm --o $output --regheader --no-save-reg


# binarize MNI 2mm mask
for ((i=2;i<=2;i++))
do
thres=`echo ${i} | awk '{printf "%.1f", $1/10}' `
echo $thres
input=${mask_MNI1mm_sm}_MNI2mm.nii.gz
output=${mask_MNI1mm_sm}_MNI2mm_bin${thres}.nii.gz
mri_binarize --i $input --o $output --min ${thres}
done




