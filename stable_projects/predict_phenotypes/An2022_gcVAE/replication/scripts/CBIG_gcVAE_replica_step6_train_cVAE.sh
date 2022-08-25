#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
optim_params="hord_params.csv"
data_path=$ROOTDIR"/data/"$2"/harm_input/"$1
checkpoint_path=$ROOTDIR"/checkpoints/"$2"/harm_model/cVAE/"$1
# get optimal parameter for each fold
declare -i row=0
params_lr=(0.0)
params_drop_out=(0.0)
params_alpha=(0.0)
params_lambda_=(0.0)
params_gamma=(0.0)
params_lrstep=(0)
params_latent_dim=(0)
params_h1=(0)
params_h2=(0)
params_h3=(0)
params_h4=(0)
params_nb_layers=(0)
while IFS=$',' read -r -a array; do
    params_lr[$row]=${array[1]}
    params_drop_out[$row]=${array[2]}
    params_alpha[$row]=${array[3]}
    params_lambda_[$row]=${array[4]}
    params_gamma[$row]=${array[5]}
    params_lrstep[$row]=${array[6]}
    params_latent_dim[$row]=${array[7]}
    params_h1[$row]=${array[8]}
    params_h2[$row]=${array[9]}
    params_h3[$row]=${array[10]}
    params_h4[$row]=${array[11]}
    params_nb_layers[$row]=${array[12]}
    row+=1
done <$checkpoint_path"/"$optim_params
params_lr=("${params_lr[@]:1}")
params_drop_out=("${params_drop_out[@]:1}")
params_alpha=("${params_alpha[@]:1}")
params_lambda_=("${params_lambda_[@]:1}")
params_gamma=("${params_gamma[@]:1}")
params_lrstep=("${params_lrstep[@]:1}")
params_latent_dim=("${params_latent_dim[@]:1}")
params_h1=("${params_h1[@]:1}")
params_h2=("${params_h2[@]:1}")
params_h3=("${params_h3[@]:1}")
params_h4=("${params_h4[@]:1}")
params_nb_layers=("${params_nb_layers[@]:1}")
# train cVAE model
for i in {0..9}; do
    python -m harmonization.cVAE.train_cVAE \
        --GPU -1 \
        --data_path $data_path/$i \
        --checkpoint_path $checkpoint_path/$i \
        --isSaving True \
        --lr ${params_lr[$i]} \
        --drop_out ${params_drop_out[$i]} \
        --alpha ${params_alpha[$i]} \
        --lambda_ ${params_lambda_[$i]} \
        --gamma ${params_gamma[$i]} \
        --lr_step ${params_lrstep[$i]} \
        --latent_dim ${params_latent_dim[$i]} \
        --nb_layers ${params_nb_layers[$i]} \
        --h1 ${params_h1[$i]} \
        --h2 ${params_h2[$i]} \
        --h3 ${params_h3[$i]} \
        --h4 ${params_h4[$i]}
done
