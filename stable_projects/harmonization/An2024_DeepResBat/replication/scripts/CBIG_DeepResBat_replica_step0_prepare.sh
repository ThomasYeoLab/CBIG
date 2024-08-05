#!/bin/bash

# Step 0
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
REPDATADIR=$CBIG_REPDATA_DIR"/stable_projects/harmonziation/An2024_DeepResBat"
cd $ROOTDIR

# create folders
mkdir -p $ROOTDIR"/job_logs"
mkdir -p $ROOTDIR"/checkpoints"
mkdir -p $ROOTDIR"/results"
mkdir -p $HOME/.cache/Python

# copy raw data and searched hyperparameters
rsync -r $REPDATADIR"/data/splits" $ROOTDIR"/data"
rsync -r $REPDATADIR"/checkpoints" $ROOTDIR
