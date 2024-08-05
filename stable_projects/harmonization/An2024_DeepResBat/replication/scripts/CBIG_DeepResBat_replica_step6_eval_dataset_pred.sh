#!/bin/bash

# Step 6
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
source activate CBIG_An2024
module load cuda/11.0
cd $ROOTDIR
export PYTHONPYCACHEPREFIX="${HOME}/.cache/Python"

# set paths
checkpoint_path=$ROOTDIR"/checkpoints/unmatch2match/eval_model/DatasetPred/"$1
output_path=$ROOTDIR"/results/dataset_prediction/"$1

if [[ "$2" = "unharm" ]]; then
    # 1. for unharmoinzed
    data_path=$ROOTDIR"/data/splits/"$1
    python -m evaluation.dataset_prediction.dataset_pred \
        --data_path $data_path \
        --checkpoint_path $checkpoint_path \
        --output_path $output_path \
        --train_name "unmatch2match_train" \
        --val_name "unmatch2match_val" \
        --test_name "unmatch2match_test" \
        --sufix "unharm"
elif [[ "$2" = "ComBat4cov" ]]; then
    # 2. for ComBat harmonization (AGE+SEX+MMSE+DX)
    data_path=$ROOTDIR"/data/unmatch2match/harm_output/ComBat/"$1"/AGE_SEX_MMSE_DX_noref"
    python -m evaluation.dataset_prediction.dataset_pred \
        --data_path $data_path \
        --checkpoint_path $checkpoint_path \
        --output_path $output_path \
        --train_name "unmatch2match_train" \
        --val_name "unmatch2match_val" \
        --test_name "unmatch2match_test" \
        --sufix "ComBat4cov"

elif [[ "$2" = "CovBat4cov" ]]; then
    # 2. for ComBat harmonization (AGE+SEX+MMSE+DX)
    data_path=$ROOTDIR"/data/unmatch2match/harm_output/CovBat/"$1"/AGE_SEX_MMSE_DX_noref"
    python -m evaluation.dataset_prediction.dataset_pred \
        --data_path $data_path \
        --checkpoint_path $checkpoint_path \
        --output_path $output_path \
        --train_name "unmatch2match_train" \
        --val_name "unmatch2match_val" \
        --test_name "unmatch2match_test" \
        --sufix "CovBat4cov"
else
    # VAE models
    data_path=$ROOTDIR"/data/unmatch2match/harm_output/"$2"/"$1
    python -m evaluation.dataset_prediction.dataset_pred \
        --data_path $data_path \
        --checkpoint_path $checkpoint_path \
        --output_path $output_path \
        --train_name "unmatch2match_train-intermediate" \
        --val_name "unmatch2match_val-intermediate" \
        --test_name "unmatch2match_test-intermediate" \
        --sufix "$2"
fi
