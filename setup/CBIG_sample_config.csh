#! /bin/csh
#

# DO NOT CHANGE: This clears old freesurfer variables if they previously exists
if( $?FREESURFER_HOME ) then
	source $FREESURFER_HOME/bin/clear_fs_env.csh 
endif

# PLEASE CHANGE: Please specify location of CBIG repository
setenv CBIG_CODE_DIR /data/users/ngohgia/storage/CBIG

# PLEASE CHANGE: define locations for these libraries
setenv FREESURFER_HOME /apps/arch/Linux_x86_64/freesurfer/5.3.0
setenv CBIG_MATLAB_DIR /apps/arch/Linux_x86_64/matlab/R2018b
setenv CBIG_SPM_DIR    /apps/arch/Linux_x86_64/spm/spm12
setenv CBIG_AFNI_DIR   /apps/arch/Linux_x86_64/afni/AFNI_2011_12_21_1014/linux_openmp_64
setenv CBIG_ANTS_DIR   /apps/arch/Linux_x86_64/ants/ants_v2.2.0/BUILD/bin/
setenv CBIG_WB_DIR     /apps/arch/Linux_x86_64/HCP/workbench-1.1.1/
setenv CBIG_FSLDIR     /apps/arch/Linux_x86_64/fsl/5.0.8

# DO NOT CHANGE: set up your environment with the configurations above
set SETUP_PATH = $CBIG_CODE_DIR/setup/CBIG_generic_setup.csh
source $SETUP_PATH

# DO NOT CHANGE: set up temporary directory for MRIread from FS6.0 for CBIG
# members using the HPC. Other users should comment this out.
setenv TMPDIR /tmpstore

# Do NOT CHANGE: set up MATLABPATH so that MATLAB can find startup.m in our repo 
setenv MATLABPATH = $CBIG_CODE_DIR/setup

# specified the default Python environment.
# Please UNCOMMENT if you follow CBIG's set up for Python environments.
# We use Python version 3.5 as default.
# Please see $CBIG_CODE_DIR/setup/python_env_setup/README.md for more details.
# source activate CBIG_py3
