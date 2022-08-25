#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
echo "Sample data from 10% - 90% to perform sample size effect analysis"
python -m utils.data_sample \
    --root_path $ROOTDIR \
    --data_path $ROOTDIR'/data' \
    --dataset_pair $1
