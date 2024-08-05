#!/bin/bash

# Example of HORD search
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
source activate CBIG_An2024
module load cuda/11.0
cd $ROOTDIR
export PYTHONPYCACHEPREFIX="${HOME}/.cache/Python"

#setup paths
data_path=$root_dir'/examples/data/unmatch2match/harm_input/ADNI-AIBL'
harm_output_path=$root_dir'/examples/data/unmatch2match/harm_output/DeepResBat/ADNI-AIBL'
checkpoint_path=$root_dir'/examples/checkpoints/unmatch2match/harm_model/DeepResBat/ADNI-AIBL'

#computation
echo $PYTHONPATH
export NUM_MKL=1
python -m HORD.hord_search_DeeResBat \
    --executable 'python -m HORD.hord_objective_DeepResBat' \
    --data_path ${data_path}/0/ \
    --harm_output_path ${harm_output_path}/0/ \
    --checkpoint_path ${checkpoint_path}/0/ \
    --ysf_sufix _G \
    --maxeval 200 \
    --nthreads 4 \
    --nGPUs 1 \
    --log_file ${checkpoint_path}/0/hord_stats.csv
