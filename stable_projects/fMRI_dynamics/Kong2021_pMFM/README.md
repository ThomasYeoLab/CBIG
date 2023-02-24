# Parametric Mean Field Model (pMFM) Analysis
# REFERENCE
********
# BACKGROUND
We developed a spatially heterogeneous large-scale dynamical circuit model that allowed for variation in local circuit properties across the human cortex. An advanced machine learning algorithm was used to estimate the model parameters from human resting-state fMRI data. We showed that parameterizing local circuit properties with both anatomical and functional gradients was necessary for generating realistic static and dynamical functional connectivity.
# CODE RELEASE
## Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: https://github.com/HeavenBluer/Parametric-MFM-Project

## Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

To download the version of the code that was last tested, you can either

* visit this link: https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.21.1-Update_Kong2021_pMFM

or

* run the following command, if you have Git installed

`git checkout -b Kong2021_pMFM v.0.21.1-Update_Kong2021_pMFM`

# USAGE
## Setup
1. Make sure you have GPU on your machine. The whole project will utilize GPU for computation. Suggested GPU runtime CUDA version: 9.1.85
2. This code is compatible with Python 3.6, to create a Python environment similar to what was used for this project:<br />
    i. Install [Anaconda](https://www.anaconda.com/products/individual#Downloads). Choose the miniconda version based on your OS.<br />
    ii. Create Anaconda environment from our `$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/replication/config/CBIG_pMFM_python_env.yml` file by `conda env create -f $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/replication/config/CBIG_pMFM_python_env.yml`. This should create an enviroment which is same as what we used in the paper.<br />
    iii. The installation time of Anaconda and environment would take around 2-3 hours to complete.<br />
3. Make sure you have installed: Matlab 2018b
4. Follow `$CBIG_CODE_DIR/setup/README.md` to setup the config file. Instead of using `$CBIG_CODE_DIR/setup/CBIG_sample_config.sh`, you need to use `$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/replication/config/CBIG_pMFM_config.sh`

## Examples
* `examples` this folder provides a toy example for our pMFM parameter estimation process
* Please refer to the README.md file in the `example` folder for the instructions to run the example

## Unit tests
* `unit_tests` this folder runs codes in `examples` and check with the reference output

## Parametric Mean Field Model Data Preparation
* `part0_pMFM_data_preparation` this folder contains scripts to generate the parcellated HCP time serise, FC and FCD
* Please follow the README.md file this folder for the detailed information of each scripts

## Parametric Mean Field Model Main Results
* `part1_pMFM_main` this folder contains scripts to generate main quantitative results in the paper.
* Please follow the README.md file this folder for the detailed information of each scripts

## Parametric Mean Field Model Control Analysis Results
* `part2_pMFM_control_analysis` this folder contains scripts to generate control analysis quantitaive results in the paper.
* Please follow the README.md file this folder for the detailed information of each scripts


# UPDATES
Release v0.21.1 (16/9/2021): Release code for revision
Release v0.19.1 (14/3/2021): Release replication of Kong2021_pMFM
Release v0.19.0 (9/3/2021): Initial release of Kong2021_pMFM



# BUGS and QUESTIONS
Please contact Kong Xiaolu at kxl920327@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
