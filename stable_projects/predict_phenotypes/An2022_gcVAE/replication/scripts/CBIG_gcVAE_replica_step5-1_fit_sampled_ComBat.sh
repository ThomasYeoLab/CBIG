#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
out_folder=$2"perc_seed"$3
data_path=$ROOTDIR"/data/sample_size/"$out_folder"/splits/"$1
output_path=$ROOTDIR"/data/sample_size/"$out_folder"/harm_output/"$1"/ComBat"

python -m harmonization.ComBat.wrapper \
    --dataset_pair $1 \
    --exp $2 \
    --data_path $data_path \
    --output_path $output_path \
    --train_file "unmatch2match_train" \
    --val_file "unmatch2match_val"
