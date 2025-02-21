#! /bin/sh

# This config is used in the process to shrink packages
# by removing redundant files not used in the fMRI preprocessing pipeline.
# In summary, the process involves:
# 1. running fMRI preprocessing unit tests with this config
#    and using strace to identify files accessed during execution.
# 2. The identified files are then moved to a temporary directory
#    to test if the setup still works without them.
# 3. If the setup still works, we have successfully identified redundant files.
#    If the setup does not work, we can use this script move the files back
#    to their original location and try again with other identified files
# For more details, please refer to the README.md in the current directory.
#
# Note that we cannot strace within containers.
# Therefore, we will use this config to set up the same environment as the containers
# and strace the preprocessing pipeline with this environment on the HPC.
#
# Specifically, compared to the original version ($CBIG_CODE_DIR/setup/CBIG_sample_config.sh), this script:\
# 1. Sets up the environment to use packages resided in the ../apps directory.
# 2. Define variables for the MATLAB Compiler runtime.
#
# Written by Tian Fang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# DO NOT CHANGE: This clears old freesurfer variables if they previously exists
if [ -n "$FREESURFER_HOME" ]; then
    # Save LD_LIBRARY_PATH to restore it later. This variable defines shared libraries used
    # by matlab to open GUI/compute/etc.
    # clear_fs_env.csh removes it when clearing old freesurfer variables but does not set it back again.
    LD_LIBRARY_PATH_CURRENT=$LD_LIBRARY_PATH

    # Clear old freesurfer variables
    while read cmd var; do
        if [[ $cmd == "unsetenv" ]]; then
            eval "unset $var"
        fi
    done <$FREESURFER_HOME/bin/clear_fs_env.csh

    # Restore old LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CURRENT
    unset LD_LIBRARY_PATH_CURRENT
fi

# PLEASE CHANGE: Please specify location of CBIG repository
export CBIG_CODE_DIR=$HOME/storage/CBIG_private

# PLEASE CHANGE: define locations for these libraries
export APPS_DIR=${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/apps
export FREESURFER_HOME=${APPS_DIR}/freesurfer/5.3.0
# ! liscened MATLAB is only used for running the unit test code, which will
# ! invoke and use the compiled version of the preprocessing pipeline using the MATLAB runtime.
export CBIG_MATLAB_DIR=/apps/matlab/R2018b
export CBIG_SPM_DIR=/apps/spm/spm12
export CBIG_AFNI_DIR=${APPS_DIR}/afni/AFNI_2011_12_21_1014/linux_openmp_64
export CBIG_ANTS_DIR=${APPS_DIR}/ants/ants_v2.2.0/BUILD/bin/
export CBIG_WB_DIR=/apps/HCP/workbench-1.1.1/
export CBIG_FSLDIR=${APPS_DIR}/fsl/5.0.10

# DO NOT CHANGE: define locations for unit tests data and replication data
export CBIG_TESTDATA_DIR=/mnt/isilon/CSC1/Yeolab/CodeMaintenance/UnitTestData
export CBIG_REPDATA_DIR=/mnt/isilon/CSC1/Yeolab/CodeMaintenance/ReplicationData

# DO NOT CHANGE: define scheduler location
export CBIG_SCHEDULER_DIR=/opt/pbs/bin

# DO NOT CHANGE: set up your environment with the configurations above
SETUP_PATH=$CBIG_CODE_DIR/setup/CBIG_generic_setup.sh
source $SETUP_PATH

# DO NOT CHANGE: set up temporary directory for MRIread from FS6.0 for CBIG
# members using the HPC, Other users should comment this out
export TMPDIR=/tmp

# Do NOT CHANGE: set up MATLABPATH so that MATLAB can find startup.m in our repo
export MATLABPATH=$CBIG_CODE_DIR/setup

# specified the default Python environment.
# Please UNCOMMENT if you follow CBIG's set up for Python environments.
# We use Python version 3.5 as default.
# Please see $CBIG_CODE_DIR/setup/python_env_setup/README.md for more details.
source activate CBIG_py3

###############################################################################
# The following section setup the MATLAB runtime environment
###############################################################################
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/runtime/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/bin/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/sys/os/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/sys/opengl/lib/glnxa64
export MATLAB_RUNTIME_DIR=/apps/matlab/MCR_R2018b/v95
