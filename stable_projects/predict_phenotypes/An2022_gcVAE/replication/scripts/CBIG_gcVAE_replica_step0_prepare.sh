#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
REPDATADIR=$CBIG_REPDATA_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
cd $ROOTDIR

# create folders
mkdir -p $ROOTDIR"/job_logs"
mkdir -p $ROOTDIR"/data"
mkdir -p $ROOTDIR"/checkpoints"
mkdir -p $ROOTDIR"/results"

# copy raw data and searched hyperparameters
rsync -r $REPDATADIR"/data/splits" $ROOTDIR"/data"
rsync -r $REPDATADIR"/data/raw_data" $ROOTDIR"/data"
rsync -r $REPDATADIR"/checkpoints" $ROOTDIR
