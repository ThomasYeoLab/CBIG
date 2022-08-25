#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
optim_params="hord_params.csv"
data_path=$ROOTDIR"/data/$2/goalDNN_input/"$1
checkpoint_path=$ROOTDIR"/checkpoints/$2/eval_model/goalDNN/"$1
# read optimal hyper-parameters
declare -i row=0
params_lr=(0.0)
params_drop_out=(0.0)
params_lambda_dx=(0.0)
params_lambda_mmse=(0.0)
params_lrstep=(0)
params_h1=(0)
params_h2=(0)
params_h3=(0)
params_h4=(0)
params_h5=(0)
params_nb_layers=(0)
while IFS=$',' read -r -a array; do
    params_lr[$row]=${array[1]}
    params_drop_out[$row]=${array[2]}
    params_lambda_dx[$row]=${array[3]}
    params_lambda_mmse[$row]=${array[4]}
    params_lrstep[$row]=${array[5]}
    params_h1[$row]=${array[6]}
    params_h2[$row]=${array[7]}
    params_h3[$row]=${array[8]}
    params_h4[$row]=${array[9]}
    params_h5[$row]=${array[10]}
    params_nb_layers[$row]=${array[11]}
    row+=1
done <$checkpoint_path"/"$optim_params
params_lr=("${params_lr[@]:1}")
params_drop_out=("${params_drop_out[@]:1}")
params_lambda_dx=("${params_lambda_dx[@]:1}")
params_lambda_mmse=("${params_lambda_mmse[@]:1}")
params_lrstep=("${params_lrstep[@]:1}")
params_h1=("${params_h1[@]:1}")
params_h2=("${params_h2[@]:1}")
params_h3=("${params_h3[@]:1}")
params_h4=("${params_h4[@]:1}")
params_h5=("${params_h5[@]:1}")
params_nb_layers=("${params_nb_layers[@]:1}")
# train goalDNN model
for i in {0..9}; do
    python -m evaluation.goalDNN.train_goalDNN \
        --GPU -1 \
        --data_path $data_path/$i \
        --checkpoint_path $checkpoint_path/$i \
        --isSaving True \
        --lr ${params_lr[$i]} \
        --drop_out ${params_drop_out[$i]} \
        --lambda_dx ${params_lambda_dx[$i]} \
        --lambda_mmse ${params_lambda_mmse[$i]} \
        --lr_step ${params_lrstep[$i]} \
        --nb_layers ${params_nb_layers[$i]} \
        --h1 ${params_h1[$i]} \
        --h2 ${params_h2[$i]} \
        --h3 ${params_h3[$i]} \
        --h4 ${params_h4[$i]} \
        --h5 ${params_h5[$i]}
done
