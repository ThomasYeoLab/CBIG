#! /bin/sh
# Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

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
    done < $FREESURFER_HOME/bin/clear_fs_env.csh

    # Restore old LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CURRENT
    unset LD_LIBRARY_PATH_CURRENT
fi

# PLEASE CHANGE: Please specify location of CBIG repository
export CBIG_CODE_DIR=$HOME/storage/CBIG_private

# PLEASE CHANGE: Please specify location of the replication input data
export CBIG_EPILEPSY_DATA_DIR=/mnt/nas/CSC11/Yeolab/Data/PNN_Lab/datasets/Lim2026_MSHBM_Epilepsy

# PLEASE CHANGE: define locations for these libraries
export FREESURFER_HOME=/apps/freesurfer/7.4.1
export CBIG_MATLAB_DIR=/apps/matlab/R2025b
export CBIG_SPM_DIR=/apps/spm/spm25
export CBIG_AFNI_DIR=/apps/afni/AFNI_25.2.18/linux_openmp_64
export CBIG_ANTS_DIR=/apps/ants/ants_v2.6.3/bin
export CBIG_WB_DIR=/apps/HCP/workbench-2.1.0
export CBIG_FSLDIR=/apps/fsl/6.0.7.19

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
source /home/mervyn/storage/Toolbox/miniconda3/etc/profile.d/conda.sh
conda activate CBIG_py3
