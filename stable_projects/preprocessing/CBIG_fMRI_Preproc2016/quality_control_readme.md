# Quality control

## T1 - T2 registration:
	
- View the intra-subject registration cost `${sub_dir}/${subject}/qc/CBIG_preproc_bbregister_intra_sub_reg.cost`. It is acceptable if the number is lower than 0.5. If it is too large, e.g. 0.7 or 0.8, this subject should be discard.
 
- Visualization check by tkregister2. The user can check the registration by looking at the alignment of the green line and the boundary of grey matter, and also by clicking "compare" to load in anatomical image and compare. The user is able to select different space coordinates to check the alignment in different position. The following three lines are the commands that the user would need. (`REG_stem` is the stem that passed into registration step, i.e. the stem of the image at which step that the user wants to use for registration)

```
cd ${sub_dir}/${subject}/bold/${run_folder}
setenv SUBJECTS_DIR ${anat_d}
tkregister2 --mov ${subject}_bld${run_folder}${REG_stem}.nii.gz --reg init-fsl/${subject}_bld${run_folder}${REG_stem}_reg.dat --surf
```
	  
----

## Motion correction and motion scrubbing (censoring)

For the parameters and interpolation method used for censoring, please refer to:
Power et al. "Methods to detect, characterize, and remove motion artifact in resting state fMRI." Neuroimage 84 (2014): 320-341.

- Motion outliers text file.

  In this file, each line is a number 0 or 1 corresponding to one frame. 0 represents this frame is censored (removed in analyses). 1 means this frame is uncensored (kept).
  `${sub_dir}/${subject}/qc/${subject}_bld${run_folder}_FDRMS0.2_DVARS50_motion_outliers.txt`

- Motion correction parameters figures.

  Please take a look at the framewise displacement (FD) and DVARS figure to see if too many frames exceed the threshold, which is:
  `${sub_dir}/${subject}/qc/${subject}_bld${run_folder}_FDRMS${FD_thres}_DVARS${DVARS_thres}.jpg`

  Another figure the user can refer to is the correlation between FD and DVARS, to see if these two parameters are consistent:
  `${sub_dir}/${subject}/qc/${subject}_bld${run_folder}_DVARS_FDRMS_correlation.jpg` (correlation between FD and DVARS)

- Motion correction quantitative numbers

  The user can check the absolute/relative mean value of motion: 
  `${sub_dir}/${subject}/qc/${subject}_bld${run_folder}_mc_abs_mean.rms`
  `${sub_dir}/${subject}/qc/${subject}_bld${run_folder}_mc_rel_mean.rms`

- Censoring plots

  In censoring QC step, the voxels are sorted by correlation or fractional difference (both in ascending ordering) between the original signal and the final censored signal. The voxels with ranking 0%, 20%, 40%, 60%, 80%, and 100% are picked up, and three signals of these voxels are plotted, which are: the original signal, the recovered signal (intermediate result), and the combined signal (replace the original signal with the recovered one at "bad" frames). The procedure are repeated both within whole brain mask and grey matter mask. The following are the figures that users can check:

  The three signals of 6 voxels ranked by correlation within grey matter mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_corr_GM_6plots.jpg`

  The three signals of 6 voxels ranked by fractional difference within grey matter mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_FracDiff_GM_6plots.jpg`

  The three signals of 6 voxels ranked by correlation within whole brain mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_corr_whole_6plots.jpg`

  The three signals of 6 voxels ranked by fractional difference within whole brain mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_FracDiff_whole_6plots.jpg`

- Censoring numerical check

  Some statistics of the correlation and fractional difference described in "Censoring plots" section. The users can check the min, max, mean, and median of the correlation and fractional difference. The text files are:

  Min, max, mean, and median correlation within grey matter mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_corr_GM.txt`

  Min, max, mean, and median fracitonal difference within grey matter mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_FracDiff_GM.txt`

  Min, max, mean, and median correlation within whole brain mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_corr_whole.txt`

  Min, max, mean, and median fractional difference within whole brain mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_FracDiff_GM.txt`
	
- Censoring correlation and fractional difference volumes

  If the users want to further know the most different regions before and after censoring, they can visualize by:

  Correlation within grey matter mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_corr_GM.nii.gz`

  Fractional difference within grey matter mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_FracDiff_GM.nii.gz`

  Correlation within whole brain mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_corr_whole.nii.gz`

  Fractional difference within whole brain mask:
  `${sub_dir}/${subject}/qc/censor_interp/${subject}_bld${run_folder}_interp_FracDiff_GM.nii.gz`
	
## Volumetric projection (inter-subject registration):

For checking purpose, our pipeline also projects the subject anatomical image to MNI 1mm space. To check if the volumetric projection (if any) has been done properly, the user can load in the projected anatomical image and the MNI 1mm template together to see if they are aligned. Here is the command:
```
freeview ${sub_dir}/${subject}/vol/norm_MNI152_1mm.nii.gz ${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz
```
