## Volumetric Templates

**MNI152**

`FSL_MNI152_FS4.5.0` contains the recon-all outputs of the 1mm MNI152 template from FSL v5.0.8. 

- This is the 6th generation template, formed by nonlinearly registering 152 individual brains and averaging across them.

- The original template can be found at `FSL_MNI152_FS4.5.0/mri/orig/001.mgz`. 
The fullly processed normalised volume can be found at `FSL_MNI152_FS4.5.0/mri/norm.mgz`.

- The reconstructed surface can be found in `FSL_MNI152_FS4.5.0/surf` 


**Colin27 (or MNI single subject)**

`SPM_Colin27_FS4.5.0` contains the recon-all outputs of the 1mm Colin27 template from SPM Anatomy Toolbox v2.2c. 

- This is a high-resolution average image across 27 scans of one individual' brain.

- The original template can be found at `SPM_Colin27_FS4.5.0/mri/orig/001.mgz`. 
The fullly processed normalised volume can be found at `SPM_Colin27_FS4.5.0/mri/norm.mgz`.

- The reconstructed surface can be found in `SPM_Colin27_FS4.5.0/surf`

## Usage

To use recon-all outputs in any FreeSurfer function:

1. set subject directory to this folder using (for csh)
```csh
setenv SUBJECTS_DIR $CBIG_CODE_DIR/data/templates/volume
```
   or (for bash)
```bash
export SUBJECTS_DIR=$CBIG_CODE_DIR/data/templates/volume
```
2. The folder name can be used as the subject ID to be passed into the command. 

   For example, to project data from MNI152 space to fsaverage (left hemisphere) by samping from the middle of the cortex
```
mri_vol2surf --mov input.nii.gz --hemi lh --projfrac 0.5 --trgsubject fsaverage --regheader FSL_MNI152_FS4.5.0 --o output.nii.gz --reshape
```
