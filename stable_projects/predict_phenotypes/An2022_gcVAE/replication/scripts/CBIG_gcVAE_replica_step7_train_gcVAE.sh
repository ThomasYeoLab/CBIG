#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
optim_params="grid_params.csv"
cVAE_model_path=$ROOTDIR"/checkpoints/$2/harm_model/cVAE/"$1
goalDNN_model_path=$ROOTDIR"/checkpoints/$2/eval_model/goalDNN/"$1
harm_input_path=$ROOTDIR"/data/$2/harm_input/"$1
goalDNN_input_path=$ROOTDIR"/data/$2/goalDNN_input/"$1
checkpoint_path=$ROOTDIR"/checkpoints/$2/harm_model/gcVAE/"$1
# get optimal parameter for each fold
declare -i row=0
params_lr=(0.0)
params_dx=(0.0)
params_mmse=(0.0)
while IFS=$',' read -r -a array; do
    params_lr[$row]=${array[1]}
    params_dx[$row]=${array[2]}
    params_mmse[$row]=${array[3]}
    row+=1
done <$checkpoint_path"/"$optim_params
params_lr=("${params_lr[@]:1}")
params_dx=("${params_dx[@]:1}")
params_mmse=("${params_mmse[@]:1}")
# train gcVAE model
for i in {0..9}; do
    python -m harmonization.gcVAE.train_gcVAE \
        --GPU -1 \
        --lr ${params_lr[$i]} \
        --lambda_dx ${params_dx[$i]} \
        --lambda_mmse ${params_mmse[$i]} \
        --cVAE_model $cVAE_model_path/$i/cVAE.pt \
        --goalDNN_model $goalDNN_model_path/$i/goalDNN.pt \
        --harm_input_path $harm_input_path/$i \
        --goalDNN_input_path $goalDNN_input_path/$i \
        --checkpoint_path $checkpoint_path/$i \
        --isSaving True \
        --exp $2
done
