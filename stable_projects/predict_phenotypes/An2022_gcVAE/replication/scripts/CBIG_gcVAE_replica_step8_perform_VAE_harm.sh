#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
raw_data_path=$ROOTDIR"/data/splits"
harm_input_path=$ROOTDIR"/data/"$2"/harm_input"
harm_output_path=$ROOTDIR"/data/"$2"/harm_output"
checkpoint_path=$ROOTDIR"/checkpoints/"$2

# cVAE harmonization
python -m harmonization.cVAE.predict_cVAE \
    --raw_data_path $raw_data_path \
    --harm_input_path $harm_input_path \
    --harm_output_path $harm_output_path \
    --checkpoint_path $checkpoint_path \
    --dataset_pair $1 \
    --exp $2
# gcVAE harmonization
python -m harmonization.gcVAE.predict_gcVAE \
    --raw_data_path $raw_data_path \
    --harm_input_path $harm_input_path \
    --harm_output_path $harm_output_path \
    --checkpoint_path $checkpoint_path \
    --dataset_pair $1 \
    --exp $2
