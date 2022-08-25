#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
out_folder=$2"perc_seed"$3
checkpoint_path=$ROOTDIR"/checkpoints/sample_size/"$out_folder"/eval_model/XGBoost/"$1
output_path=$ROOTDIR"/results/sample_size/"$out_folder"/eval_model/XGBoostPred/"$1

# 1. for ComBat harmonization (AGE+SEX)
data_path=$ROOTDIR"/data/sample_size/"$out_folder"/harm_output/"$1"/ComBat/AGE_SEX_noref"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name "unmatch2match_train_full" \
    --val_name "unmatch2match_val_full" \
    --test_name "unmatch2match_test" \
    --save_suffix "ComBat"
# 2. for cVAE harmonization
data_path=$ROOTDIR"/data/sample_size/"$out_folder"/harm_output/"$1"/cVAE"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name "unmatch2match_train_full-intermediate" \
    --val_name "unmatch2match_val_full-intermediate" \
    --test_name "unmatch2match_test-intermediate" \
    --save_suffix "cVAE"
# 3. for gcVAE harmonization
data_path=$ROOTDIR"/data/sample_size/"$out_folder"/harm_output/"$1"/gcVAE"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name "unmatch2match_train_full-intermediate" \
    --val_name "unmatch2match_val_full-intermediate" \
    --test_name "unmatch2match_test-intermediate" \
    --save_suffix "gcVAE"
