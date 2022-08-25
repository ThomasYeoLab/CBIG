#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# 1. generated input file for goalDNN
out_folder=$2"perc_seed"$3
data_path=$ROOTDIR"/data/sample_size/"$out_folder
python -m utils.input_generation \
    --data_path $data_path \
    --dataset_pair $1 \
    --exp "sample_size" \
    --harm True
# 2. goalDNN makes prediction
data_path=$ROOTDIR"/data/sample_size/"$out_folder"/goalDNN_input/"$1
# using goalDNN devloped on full unmatch2match set to evaluate
checkpoint_path=$ROOTDIR"/checkpoints/unmatch2match/eval_model/goalDNN/"$1
save_path=$ROOTDIR"/results/sample_size/"$out_folder"/eval_model/goalDNNPred/"$1
python -m evaluation.goalDNN.predict_goalDNN \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --save_path $save_path \
    --dataset_pair $1 \
    --exp "sample_size"
# 3. evaluate prediction
python -m evaluation.goalDNN.eval_prediction \
    --save_path $save_path \
    --dataset_pair $1 \
    --exp "sample_size"
