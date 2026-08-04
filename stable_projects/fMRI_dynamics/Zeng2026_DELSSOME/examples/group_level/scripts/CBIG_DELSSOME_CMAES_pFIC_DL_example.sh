#!/bin/sh
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

source activate CBIG_DELSSOME

example_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
python_dir=${example_dir}/../../group_level

# select mode
mode=train_val_test

# path to the config file (includes pretrained predictor model path)
config_path=${example_dir}/config/example_CMAES_pFIC_DELSSOME.ini

# number of random seeds used
n_seed=1

# path to the folder with all the input files
input_dir=${example_dir}/input

# path to the folder under which all the output files will be saved
output_dir=${example_dir}/output/dl

# specify the total number of training epoch for parameter optimization
num_epochs=5

# the number of trials to select new seeds
trials=1

# the training mode. DELSSOME scores each candidate parameter set with the deep
# learning cost predictor instead of Euler integration; Euler does the integration
train_mode=DELSSOME

# run the example
cd ${python_dir}
python -u CBIG_DELSSOME_CMAES_main_pFIC.py --mode ${mode} --config_path ${config_path} --n_seed ${n_seed} \
    --input_dir ${input_dir} --output_dir ${output_dir} --num_epochs ${num_epochs} \
    --trials ${trials} --train_mode ${train_mode} --seed 123456

# check example output
cd ${example_dir}
python -u CBIG_DELSSOME_check_CMAES_pFIC_DL_example_results.py

conda deactivate
