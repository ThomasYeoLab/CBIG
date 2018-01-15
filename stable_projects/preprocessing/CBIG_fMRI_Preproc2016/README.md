## Reference

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

----

## Code Release

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
  
  The compulsory softwares include FreeSurfer (5.3 or 4.5), FSL (5.0.8), Matlab (2014a), and Python (2 or 3, only build-in functions are needed). If the user wants to use `CBIG_preproc_despiking` step, then AFNI is needed. If the user wants to use `CBIG_preproc_native2mni_ants` step, then ANTs (2.2.0) is needed.

  NOTE: There is a bug in early builds of ANTs (before Aug 2014) that causes resampling for timeseries to be wrong. We have tested that our
codes would work on ANTs version 2.2.0. 

----

## Bugs and Questions

Please contact Ru(by) Kong at roo.cone@gmail.com, Jingwei Li at jingwei.li@u.nus.edu, Nanbo Sun at sun464879934@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
