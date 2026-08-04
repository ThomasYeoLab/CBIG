#!/bin/sh
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# module load cuda/12.4  # For GPU users please uncomment this line
source activate CBIG_DELSSOME

example_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
model_trainer_dir=${example_dir}/../../group_level

# path to the config file
config_path=${example_dir}/config/example_predictor_train.ini

# random seed
seed=123456

cd ${model_trainer_dir}
python -u CBIG_DELSSOME_predictor_trainer.py --config_path ${config_path} --seed ${seed}

# check example output
cd ${example_dir}
python -u CBIG_DELSSOME_check_predictor_train_results.py

conda deactivate
