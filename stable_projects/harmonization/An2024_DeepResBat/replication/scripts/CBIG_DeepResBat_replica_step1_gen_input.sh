#!/bin/bash

# Step 1
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
source activate CBIG_An2024
module load cuda/11.0
cd $ROOTDIR
export PYTHONPYCACHEPREFIX="${HOME}/.cache/Python"

# run
python -m utils.input_generation \
    --dataset_pair $1 \
    --input_path $ROOTDIR'/data/splits' \
    --harm_input_path $ROOTDIR'/data/unmatch2match/harm_input'
