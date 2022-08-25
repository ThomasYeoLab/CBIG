#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# 1. generated input file for goalDNN
data_path=$ROOTDIR"/data"
python -m utils.input_generation \
    --data_path $data_path \
    --dataset_pair $1 \
    --exp $2 \
    --harm True
# 2. goalDNN makes prediction
data_path=$ROOTDIR"/data/$2/goalDNN_input/"$1
checkpoint_path=$ROOTDIR"/checkpoints/$2/eval_model/goalDNN/"$1
save_path=$ROOTDIR"/results/"$2"/eval_model/goalDNNPred/"$1
python -m evaluation.goalDNN.predict_goalDNN \
    --data_path $data_path \
    --checkpoint_path $checkpoint_path \
    --save_path $save_path \
    --dataset_pair $1 \
    --exp $2
# 3. evaluate prediction
python -m evaluation.goalDNN.eval_prediction \
    --save_path $save_path \
    --dataset_pair $1 \
    --exp $2
