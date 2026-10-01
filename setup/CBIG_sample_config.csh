#! /bin/csh
#

# DO NOT CHANGE: This clears old freesurfer variables if they previously exists
if( $?FREESURFER_HOME ) then
    # Save LD_LIBRARY_PATH to restore it later. This variable defines shared libraries used
    # by matlab to open GUI/compute/etc.
    # clear_fs_env.csh removes it when clearing old freesurfer variables but does not set it back again.
    set LD_LIBRARY_PATH_CURRENT=$LD_LIBRARY_PATH

    # Clear old freesurfer variables
    source $FREESURFER_HOME/bin/clear_fs_env.csh 

    # Restore old LD_LIBRARY_PATH
    setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH_CURRENT
    unset LD_LIBRARY_PATH_CURRENT
endif

# PLEASE CHANGE: Please specify location of CBIG repository
setenv CBIG_CODE_DIR $HOME/storage/CBIG_private

# PLEASE CHANGE: define locations for these libraries
setenv FREESURFER_HOME    /apps/freesurfer/7.4.1
setenv CBIG_MATLAB_DIR    /apps/matlab/R2025b
setenv CBIG_SPM_DIR       /apps/spm/spm25
setenv CBIG_AFNI_DIR      /apps/afni/AFNI_25.2.18/linux_openmp_64
setenv CBIG_ANTS_DIR      /apps/ants/ants_v2.6.3/bin
setenv CBIG_WB_DIR        /apps/HCP/workbench-2.1.0
setenv CBIG_FSLDIR        /apps/fsl/6.0.7.19

# DO NOT CHANGE: define locations for unit tests data and replication data
setenv CBIG_TESTDATA_DIR  /mnt/nas/CSC21/Yeolab/CodeMaintenance/UnitTestData
setenv CBIG_REPDATA_DIR   /mnt/nas/CSC21/Yeolab/CodeMaintenance/ReplicationData

# DO NOT CHANGE: define scheduler location
setenv CBIG_SCHEDULER_DIR /opt/pbs/bin

# DO NOT CHANGE: set up your environment with the configurations above
set SETUP_PATH = $CBIG_CODE_DIR/setup/CBIG_generic_setup.csh
source $SETUP_PATH

# DO NOT CHANGE: set up temporary directory for MRIread from FS6.0 for CBIG
# members using the HPC. Other users should comment this out.
setenv TMPDIR /tmp

# Do NOT CHANGE: set up MATLABPATH so that MATLAB can find startup.m in our repo 
setenv MATLABPATH $CBIG_CODE_DIR/setup

# specified the default Python environment.
# Please UNCOMMENT if you follow CBIG's set up for Python environments.
# We use Python version 3.5 as default.
# Please see $CBIG_CODE_DIR/setup/python_env_setup/README.md for more details.
# source CBIG_init_conda
# conda activate CBIG_py3
