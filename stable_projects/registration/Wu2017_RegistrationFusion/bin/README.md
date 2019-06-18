## Stand-alone Usage

To use volume-to-fsaverage mappings without the trouble of downloading our entire repository, just download the `standalone_scripts_for_MNI_fsaverage_projection.zip` file for MNI152-to-fsaverage projections or the `standalone_scripts_for_Colin_fsaverage_projection.zip` file for Colin27-to-fsaverage projections; then follow the README.md inside.

To convert between MNI152 coordinates and fsavearge vertices without the trouble of downloading our entire repository, just download the `standalone_scripts_for_MNI_fsaverage_coordinates_conversion.zip` file; then follow the README.md inside.

Note that FreeSurfer and Matlab need to be installed before the stand-alone scripts can be run.

----

## Folder Contents

This directory contains the following folders and files:

- `final_warps_FS5.3`: final mappings for RF approaches (RF-M3Z and RF-ANTs) using FreeSurfer 5.3

- `liberal_cortex_masks_FS5.3`: liberal cortical masks used in surface-to-volume projections

- `scripts_final_proj`: scripts to directly use final mappings for RF approaches (RF-M3Z and RF-ANTs). Refer to README.md in this folder for how to use these scripts.

- `standalone_scripts_for_MNI_fsaverage_coordinates_conversion` and `standalone_scripts_for_MNI_fsaverage_coordinates_conversion.zip`: stand-alone scripts for converting between MNI152 coordinates and fsavearge vertices downloading the other folders in the CBIG repository. Refer to README.md in this folder further for how to use these scripts.

    - `standalone_scripts_for_MNI_fsaverage_coordinates_conversion.zip` contains the following folders: `standalone_scripts_for_MNI_fsaverage_coordinates_conversion`, `liberal_cortex_masks_FS5.3`, and `final_warps_FS5.3` 

- `standalone_scripts_for_MNI_fsaverage_projection` and `standalone_scripts_for_MNI_fsaverage_projection.zip`: stand-alone scripts for using the final MNI152-to-fsaverage mapping without downloading the other folders in the CBIG repository. Refer to README.md in this folder further for how to use these scripts.

    - `standalone_scripts_for_MNI_fsaverage_projection.zip` contains the following folders: `standalone_scripts_for_MNI_fsaverage_projection` and `final_warps_FS5.3` 

- `standalone_scripts_for_Colin_fsaverage_projection` and `standalone_scripts_for_Colin_fsaverage_projection.zip`: stand-alone scripts for using the final Colin27-to-fsaverage mapping without downloading the other folders in the CBIG repository. Refer to README.md in this folder further for how to use these scripts.

    - `standalone_scripts_for_Colin_fsaverage_projection.zip` contains the following folders: `standalone_scripts_for_Colin_fsaverage_projection` and `final_warps_FS5.3`

- `GSP_subjectid.csv`: list of subject IDs for GSP subjects, used by RF approaches (RF-M3Z and RF-ANTs) by default

- `colin27.register.dat`: linear registration warp between Colin27 template and fsaverage volume (used by Affine approach). This is generated using FSL, similar to FreeSurfer's own mni152.register.dat.

- `images`: image files used in some READMEs

