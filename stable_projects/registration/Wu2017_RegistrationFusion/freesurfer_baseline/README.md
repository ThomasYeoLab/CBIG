## Affine Commands

Use the following commands to implement Affine approach for volume-to-surface or surface-to-volume mapping.

**Each line in my_input_list should contain one file name. Note that inputs should be in the format of .nii or .nii.gz**

- MNI152-to-fsaverage: `./CBIG_RF_FS_proj_vol2fsaverage.sh -l my_input_list -f Affine`

- Colin27-to-fsaverage: `./CBIG_RF_FS_proj_vol2fsaverage.sh -l my_input_list -f Affine -s 'SPM_Colin27_FS4.5.0'`

**Each line in my_input_list should contain one file name. Note that the inputs should be in .mat format, containing two variables `lh_label` and `rh_label`, corresponding to the input data in left and right hemispheres respectively.**

- fsaverage-to-MNI152: `./CBIG_RF_FS_proj_fsaverage2vol.sh -l my_input_list -f Affine`

- fsaverage-to-Colin27: `PARENT_PATH=$(dirname "$(readlink -f "$(pwd)")")` 
then `./CBIG_RF_FS_proj_fsaverage2vol.sh -l my_input_list -f Affine -s 'SPM_Colin27_FS4.5.0' -c $PARENT_PATH/bin/liberal_cortex_masks_FS5.3/SPM_Colin27_FS4.5.0_cortex_estimate.nii.gz`

----

## MNIsurf/Colinsurf Commands

Use the following commands to implement MNIsurf/Colinsurf approach for volume-to-surface or surface-to-volume mapping. 

**Each line in my_input_list should contain one file name. Note that inputs should be in the format of .nii or .nii.gz**

- MNI152-to-fsaverage: `./CBIG_RF_FS_proj_vol2fsaverage.sh -l my_input_list`

- Colin27-to-fsaverage: `./CBIG_RF_FS_proj_vol2fsaverage.sh -l my_input_list -s 'SPM_Colin27_FS4.5.0'`

**Each line in my_input_list should contain one file name. Note that the inputs should be in .mat format, containing two variables `lh_label` and `rh_label`, corresponding to the input data in left and right hemispheres respectively.**

- fsaverage-to-MNI152: `./CBIG_RF_FS_proj_fsaverage2vol.sh -l my_input_list`

- fsaverage-to-Colin27: `PARENT_PATH=$(dirname "$(readlink -f "$(pwd)")")` 
then `./CBIG_RF_FS_proj_fsaverage2vol.sh -l my_input_list -s 'SPM_Colin27_FS4.5.0' -c $PARENT_PATH/bin/liberal_cortex_masks_FS5.3/SPM_Colin27_FS4.5.0_cortex_estimate.nii.gz`

