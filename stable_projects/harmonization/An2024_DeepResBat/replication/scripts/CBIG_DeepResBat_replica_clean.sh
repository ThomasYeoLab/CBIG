#!/bin/bash

# Clean all
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
cd $ROOTDIR

# run
rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/data/unmatch2match

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/data/splits

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/checkpoints

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/results

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/job_logs

find "$ROOTDIR" -type d -name '__pycache__' -exec rm -rf {} +
