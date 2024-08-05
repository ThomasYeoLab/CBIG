#!/bin/bash

# Step 8
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
source activate CBIG_An2024
module load cuda/11.0
cd $ROOTDIR
export PYTHONPYCACHEPREFIX="${HOME}/.cache/Python"

# set paths
output_path=$ROOTDIR"/results/assoc_glm/"

if [[ "$2" = "unharm" ]]; then
    # 1. for unharmoinzed
    input_path=$ROOTDIR"/data/splits/"$1
    test_name="unmatch2match_test"
elif [[ "$2" = "ComBat4cov" ]]; then
    # 2. for ComBat harmonization (AGE+SEX+MMSE+DX)
    input_path=$ROOTDIR"/data/unmatch2match/harm_output/ComBat/"$1"/AGE_SEX_MMSE_DX_noref"
    test_name="unmatch2match_test"
elif [[ "$2" = "CovBat4cov" ]]; then
    # 3. for CovBat harmonization (AGE+SEX+MMSE+DX)
    input_path=$ROOTDIR"/data/unmatch2match/harm_output/CovBat/"$1"/AGE_SEX_MMSE_DX_noref"
    test_name="unmatch2match_test"
else
    # VAE models
    input_path=$ROOTDIR"/data/unmatch2match/harm_output/"$2"/"$1
    test_name="unmatch2match_test-intermediate"
fi

python -m evaluation.association_analysis.glm \
    --input_path $input_path \
    --output_path $output_path \
    --dataset_pair $1 \
    --test_name $test_name \
    --model $2 \
    --type 'DX'
python -m evaluation.association_analysis.glm \
    --input_path $input_path \
    --output_path $output_path \
    --dataset_pair $1 \
    --test_name $test_name \
    --model $2 \
    --type 'MMSE'
