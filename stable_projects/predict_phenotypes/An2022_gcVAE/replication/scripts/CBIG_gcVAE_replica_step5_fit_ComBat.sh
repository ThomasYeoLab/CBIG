#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

data_path=$ROOTDIR"/data/splits/"$1
output_path=$ROOTDIR"/data/"$2"/harm_output/"$1"/ComBat"

python -m harmonization.ComBat.wrapper \
    --dataset_pair $1 \
    --exp $2 \
    --data_path $data_path \
    --output_path $output_path \
    --train_file $2"_train" \
    --val_file $2"_val"
