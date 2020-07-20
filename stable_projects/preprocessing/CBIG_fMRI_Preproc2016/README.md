## Reference

Li J, Kong R, Liegeois R, Orban C, Tan Y, Sun N, Holmes AJ, Sabuncu MR, Ge T, Yeo BTT, [**Global signal regression strengthens association between resting-state functional connectivity and behavior**](https://doi.org/10.1016/j.neuroimage.2019.04.016), Neuroimage, 2019, 196:126-141.

Kong R, Li J, Orban C, Sabuncu MR, Liu H, Schaefer A, Sun N, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT, [**Spatial Topography of Individual-Specific Cortical Networks Predicts Human Cognition, Personality and Emotion**](https://academic.oup.com/cercor/advance-article/doi/10.1093/cercor/bhy123/5033556?guestAccessKey=2fa23bc8-59c7-4ff1-9360-1846d472c6dd), Cerebral Cortex, 2019, 26(6):2533-2551

----

## Background

This folder contains a resting-state fMRI preprocessing pipeline written by CBIG group. Our preprocessing pipeline allows flexible preprocessing order by specifying the order of preprocessing steps in a configuration text file. The preprocessing steps include:
- slice-time correction
- motion correction
- spatial distortion correction
- intra-subject registration between T1 and T2* images
- nuisance regression
- temporal interpolation of censored frames
- bandpass filtering
- projections to standard surface & volumetric spaces
- functional connectivity (FC) matrix computation

We also provide multiple types of quality control (QC) figures for data inspection. For example:

- grey-scale timeseries intensity plot for voxels within grey matter (e.g. left subfigure)
- group-level dependency between QC-FC correlation and ROI-to-ROI Euclidean distance (e.g. right subfigure)

<img src="readme_figures/README_figure.png" height="400" />

----

## Code Release

### Download stand-alone repository

Since the whole GitHub repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link:
[https://github.com/ThomasYeoLab/Standalone_CBIG_fMRI_Preproc2016](https://github.com/ThomasYeoLab/Standalone_CBIG_fMRI_Preproc2016)


### Download whole repository

Except for this project, if you want to use the code for other stable projects from our lab as well, you need to download the whole repository.

- To download the version of the code that was last tested, you can either

  - visit this link:
  [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.17.2-CBIG_preproc_spatial_distortion_correction](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.17.2-CBIG_preproc_spatial_distortion_correction)
  
  or
  
  - run the following command, if you have Git installed
  
  ```
  git checkout -b CBIG_fMRI_Preprocessing v0.17.2-CBIG_preproc_spatial_distortion_correction
  ```

### Usage 

- How to run our preprocessing pipeline?

  The main function is `CBIG_preproc_fMRI_preprocess.csh`. Type `CBIG_preproc_fMRI_preprocess.csh -help` in command line to see the details.
  
  To use the script, the user needs to pass in a configuration file. An example of the configuration file is given in this file: `example_config.txt`.
  
  Each uncommented line (not started by #) in the configuration file corresponds to one script in this folder `$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016`. For example,
  `CBIG_preproc_skip -skip 4` corresponds to `CBIG_preproc_skip.csh`. 
  
  By using different configuration file, the user can change the ordering of preprocessing steps and the parameters passed in each step. To get more details of each step, the user can type in `./CBIG_preproc_xxx.csh -help` in command line.
  
- The utilities folder contains all matlab functions that are called by the c shell scripts. The user can type in `help function_name` in a Matlab command window to see how each function works.

- Check errors

  It is possible that some subjects are not preprocessed successfully (e.g. matlab crash/license issue). Therefore, we provide a check error function `CBIG_preproc_check_error.csh` which goes through each log file to search for the error and successful messages. Type `./CBIG_preproc_check_error.csh -help` for more information. If a subject fails. we recommend the users to delete the whole subject folder and rerun the pipeline from beginning.
  
- Quality control

  To see how to perform quality control for your preprocessed data, please check `quality_control_readme.md`.
  
- Unit test

  We provide the instruction of how to do the unit test: `unit_test_clustering_100_subjects_readme.md` in this folder. 
  Remark: Since the unit test includes data and scripts which are specific for our server, it is only applicable within CBIG lab.
  
- Software versions
  
  The compulsory softwares include FreeSurfer (5.3 or 4.5), FSL (5.0.10), Matlab (2014a), and Python (2 or 3, only build-in functions are needed). If the user wants to use `CBIG_preproc_despiking` step, then AFNI is needed. If the user wants to use `CBIG_preproc_native2mni_ants` step, then ANTs (2.2.0) is needed.

  NOTE: There is a bug in early builds of ANTs (before Aug 2014) that causes resampling for timeseries to be wrong. We have tested that our
codes would work on ANTs version 2.2.0. 

----

## Updates

- Release v0.4.4 (20/10/2017): Initial release of CBIG fMRI preprocessing pipeline.
- Release v0.4.5 (01/12/2017):

	1. Change motion correction (mcflirt) interpolation method from default **trilinear** to **spline**.
	
	2. Add an optional preprocessing step to perform despiking by **AFNI 3dDespike**.
	
	3. Add a preprocessing step to generate ROIs2ROIs functional connectivity matrix for input subject. 
	
- Release v0.4.6 (04/01/2018): 
    
    1. Add functionality: projecting fMRI data from subject-specific space to MNI 2mm space using ANTs. 
    
    2. Speed up censoring interpolation step by applying loose whole brain mask and optimizing the number of voxels processed each time. 
    
    3. Add some functionality to generate more QC plots (plots of mcflirt parameters; grey plots reflecting signal intensity in grey matter). 
    
    4. Force medial wall vertices to be NaN for the data in fsaverage surface space.
    
- Release v0.4.7 (16/01/2018):

    1. Add `config` and `unit_tests` folders.
    
    2. Add project-specific prefix `preproc` for preprocessing scripts.

- Release v0.4.9 (31/01/2018):

    1. Fix broken symbolic link under `bin` folder.
    
    2. Add preprocessing scripts to plot QC-FC correlation versus ROIs to ROIs distance.
    
    3. Add an option to specify maximal memory usage in censoring interpolation step.
    
    4. Add `examples` folder.
    
- Release v0.4.11 (21/03/2018):

    1. Fix a bug: ventricles mask was not generated when the ventricles segmentation in anantomical space is <=100 voxels.
    
    2. Change the constraints of minimal number of voxels in ventricles/wm masks to be constraints of total volume in ventricles/wm masks.
    
    3. Modify functions `CBIG_ComputeROIs2ROIsCorrelationMatrix.m` (used by `CBIG_preproc_FCmetrics.m`) and `CBIG_preproc_QC_greyplot.m` so that they do not need matlab statistics toolbox now.
    
    4. In regression step, remove the exact zero columns of regressors.
    
- Release v0.6.2 (15/07/2018):

    1. Update README.md for creating stand-alone repo.
    
    2. Add unit test: check correctness of 419x419 function connectivity matrix.
    
    3. Add scripts for all unit tests.
    
- Release v0.9.6 (12/04/2019):
    1. Update censoring interpolation step to avoid crashing when there are vetices whose timeseries are all NaNs and to reduce run time if no frame needs to be censored.
    
    2. Update QC greyplot step: global signal to be plotted is now computed from the input volume of the nuisance regression step (e.g. after T1-T2* registration and before censoring interpolation, if you use the default config file), whereas previously it was computed from the input volume of the QC greyplot step (e.g. after bandpass filtering, if you use the default config file).
    
    3. Fix a bug to avoid crashing when -final_mask option is not passed into `CBIG_preproc_native2mni.csh` and `CBIG_preproc_native2mni_ants.csh`
    
    4. Save out motion correction transformation matrices for future distortion correction usage (distortion correction scripts have not been added in this release. They will be added in a future release.)
    
    5. Change the rule-of-thumb of choosing bbregister transformation matrix. In previous versions, the registration matrix for each run was replaced with the registration matrix of the run who had the lowest BBR cost (i.e. the best run) in the same subject. From this release onwards, the BBR registration of the best run is applied to other runs only if the registration can improve the BBR cost of other runs. 
    
    6. Use the run with lowest BBR cost, instead of the first run, to create all the masks in the functional native volumetric space.
    
    7. Remove `-censor` option in `CBIG_preproc_bandpass_fft.csh` and options of `-low_f` and `-high_f` in `CBIG_preproc_censor.csh`. Include a readme about bandpass filtering and censoring (i.e. `$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/Recommendation_of_bandpass_censoring.md`).
    
- Release v0.9.8 (30/04/2019):
    1. Bug fix: the script `$CBIG_CODE_DIR/utilities/scripts/CBIG_antsReg_vol2vol.sh` was supposed to be released in v0.9.6 but not released, causing crash of `CBIG_preproc_native2mni_ants` step. In this version, the updated `$CBIG_CODE_DIR/utilities/scripts/CBIG_antsReg_vol2vol.sh` will be released.
    
- Release v0.13.1 (19/07/2019): Update references in the readme.
    
- Release v0.17.0 (19/02/2020): Avoid using absolute paths. Add new environment variables to avoid possible problems caused by hard-coded absolute paths.

- Release v0.17.2 (07/07/2020):
	1. Bug fix: Fix 'out-of-bound' error of `CBIG_preproc_fslmcflirt_outlier.csh` due to incorrect extraction of number of frames from `$boldfile"_mc_tmp.nii.gz"`.
	2. Bug fix: Fix the threshold for ventricle mask erosion in functional space from `$num_vent` (number of voxels) to `$vent_vol` (volume).
	3. Add spatial distortion correction: CBIG_preproc_spatial_distortion_correction.csh.
	4. Spatial distortion correction step requires a newer FSL version. Update the default FSL verision to 5.0.10.

----

## Bugs and Questions

Please contact Ru(by) Kong at roo.cone@gmail.com, Jingwei Li at jingweili.sjtu.nus@gmail.com, Nanbo Sun at sun464879934@gmail.com, Shaoshi Zhang at 0zhangshaoshi0@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
