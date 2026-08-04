#! /bin/sh
# Last successfully run on July 4th, 2026 with git repository version v0.1.0-Zeng2026_DELSSOME
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# DO NOT CHANGE: This clears old freesurfer variables if they previously exists
if [ -n "$FREESURFER_HOME" ]; then
    $FREESURFER_HOME/bin/clear_fs_env.csh
fi

# PLEASE CHANGE: Please specify location of CBIG repository
export CBIG_CODE_DIR=$HOME/storage/CBIG_private

# DO NOT CHANGE: define locations for unit tests data and replication data
export CBIG_TESTDATA_DIR=/mnt/isilon/CSC1/Yeolab/CodeMaintenance/UnitTestData
export CBIG_REPDATA_DIR=/mnt/isilon/CSC1/Yeolab/CodeMaintenance/ReplicationData

# DO NOT CHANGE: define scheduler location
export CBIG_SCHEDULER_DIR=/opt/pbs/bin

# DO NOT CHANGE: set up your environment with the configurations above
SETUP_PATH=$CBIG_CODE_DIR/setup/CBIG_generic_setup.sh
source $SETUP_PATH

# DO NOT CHANGE: set up temporary directory for MRIread from FS6.0 for CBIG
# members using the HPC. Other users should comment this out
export TMPDIR=/tmp

# The Zeng2026_DELSSOME project is Python-only. Create the environment once with
#   conda env create -f replication/config/CBIG_DELSSOME_python_env.yml
# and activate it before running the group-level pipeline. The individual-level
# pipeline uses a separate environment (indiv_level/environment.yml, adds R).
# source activate CBIG_DELSSOME
