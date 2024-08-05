#!/bin/bash

# Step 2
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
source activate CBIG_An2024
module load cuda/11.0
cd $ROOTDIR
export PYTHONPYCACHEPREFIX="${HOME}/.cache/Python"

# setup paths
data_path=$ROOTDIR"/data/splits/"$1

# ComBat
output_path=$ROOTDIR"/data/unmatch2match/harm_output/ComBat/"$1
python -m harmonization.ComBat.wrapper \
    --dataset_pair $1 \
    --data_path $data_path \
    --output_path $output_path
# CovBat
output_path=$ROOTDIR"/data/unmatch2match/harm_output/CovBat/"$1
python -m harmonization.ComBat.wrapper_CovBat \
    --dataset_pair $1 \
    --data_path $data_path \
    --output_path $output_path
