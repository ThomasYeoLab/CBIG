#! /bin/csh -f 
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

foreach type (Yeo2011_7Networks_N1000.split_components Yeo2011_17Networks_N1000.split_components)
    mri_vol2vol --targ $CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz --mov ../MNI152/$type.FSL_MNI152_FreeSurferConformed_1mm.nii.gz --o ../MNI152/$type.FSL_MNI152_2mm.nii.gz --regheader --no-save-reg
    mri_vol2vol --targ $CBIG_CODE_DIR/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_1mm_brain.nii.gz --mov ../MNI152/$type.FSL_MNI152_FreeSurferConformed_1mm.nii.gz --o ../MNI152/$type.FSL_MNI152_1mm.nii.gz --regheader --no-save-reg 
end
