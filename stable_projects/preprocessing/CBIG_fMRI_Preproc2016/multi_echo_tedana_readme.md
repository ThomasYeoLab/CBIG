# Multi Echo Preprocessing (ME) README

## Overview
Multi-echo (ME) fMRI refers to collecting data at multiple echo times, resulting in multiple volumes with varying levels of contrast acquired per radio frequency pulse. By collecting multi-echo data, we can compare results across different echoes, more importantly combine the results by weighted averaging and denoise the data based on information contained in the echoes. 'Tedana' is an ICA-based denoising pipeline built especially for multi-echo data. Rather than analyzing single-echo time series separately, tedana combines them into an “optimally combined time series”, and then denoising the data based on information contained in the echoes using a multi-echo ICA-based denoising method. For more detailed information about 'tedana', please refer to their website: https://tedana.readthedocs.io/en/stable/index.html


## Preparation

In order to run tedana, please make sure that the following packages are available under conda environment with name tedana:

* nilearn
* nibabel
* numpy
* scikit-learn
* scipy
* tedana
* duecredit (optional) 'duecredit' is a python package that is used, but not required by 'tedana'. These warnings do not affect any of the processing within the 'tedana'. To avoid this warning, you can install 'duecredit' with 'pip install duecredit'.

If the user has installed CBIG_py3 environment, the above list of packages are already installed. For more detailed information about our CBIG_py3 envrionment, please refer to:
https://github.com/ThomasYeoLab/CBIG/tree/master/setup/python_env_setup#quick-installation-for-linux

And also please make sure that echo times for each echo are available.


## How to run the ME script?
An example command is shown as the following:

`./CBIG_preproc_multiecho_denoise.csh -s sub005 -d ~/fmri_preprocess -bld 001 -BOLD_stem _rest_skip4_stc_mc_sdc -echo_number 3 -echo_time 12,30.11,48.22`

**[IMPORTANT]**: Please note that echo time must be in milleseconds, and echo time should be separated by comma and should be set in an ascending order.

## Quality Control (QC)

Two grey plots, which are before and after denoising, will be generated in QC of tedana. By right, the one which is after denoising should be more clean and consistent. After MEICA, global artefact should become more pronounced while local artefacts are removed. QC result will be in `$sub_dir/$subject/qc` folder.
