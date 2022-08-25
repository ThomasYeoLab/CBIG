#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
checkpoint_path=$ROOTDIR"/checkpoints/"$2"/eval_model/XGBoost/"$1
output_path=$ROOTDIR"/results/"$2"/eval_model/XGBoostPred/"$1
# 1. for unharmoinzed
data_path=$ROOTDIR"/data/splits/"$1
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name $2"_train" \
    --val_name $2"_val" \
    --test_name $2"_test" \
    --save_suffix "unharm"
# 2. for ComBat harmonization (AGE+SEX)
data_path=$ROOTDIR"/data/"$2"/harm_output/"$1"/ComBat/AGE_SEX_noref"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name $2"_train" \
    --val_name $2"_val" \
    --test_name $2"_test" \
    --save_suffix "ComBat"
# 3. for ComBat harmonization (AGE+SEX+MMSE+DX)
data_path=$ROOTDIR"/data/"$2"/harm_output/"$1"/ComBat/AGE_SEX_MMSE_DX_noref"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name $2"_train" \
    --val_name $2"_val" \
    --test_name $2"_test" \
    --save_suffix "ComBat4cov"
# 4. for cVAE harmonization
data_path=$ROOTDIR"/data/"$2"/harm_output/"$1"/cVAE"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name $2"_train-intermediate" \
    --val_name $2"_val-intermediate" \
    --test_name $2"_test-intermediate" \
    --save_suffix "cVAE"
# 5. for gcVAE
data_path=$ROOTDIR"/data/"$2"/harm_output/"$1"/gcVAE"
python -m evaluation.XGBoost.train_XGBoost \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --output_path $output_path \
    --train_name $2"_train-intermediate" \
    --val_name $2"_val-intermediate" \
    --test_name $2"_test-intermediate" \
    --save_suffix "gcVAE"
