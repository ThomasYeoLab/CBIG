#!/bin/sh
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

source activate CBIG_DELSSOME

example_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
python_dir=${example_dir}/../../group_level

# select mode
mode=train_val_test

# path to the config file
config_path=${example_dir}/config/example_CMAES_pFIC.ini

# number of random seeds used, set to a value that is > 1 will SIGNIFICANTLY slow down the script
n_seed=1

# path to the folder with all the input files
input_dir=${example_dir}/input

# path to the folder under which all the output files will be saved
output_dir=${example_dir}/output/pfic

# specify the total number of training epoch for parameter optimization
num_epochs=5

# the number of trials to select new seeds. If seed is specified, force to be 1.
trials=1

# the training mode. Choose from DELSSOME or Euler
train_mode=Euler

# run the example
cd ${python_dir}
python -u CBIG_DELSSOME_CMAES_main_pFIC.py --mode ${mode} --config_path ${config_path} --n_seed ${n_seed} \
    --input_dir ${input_dir} --output_dir ${output_dir} --num_epochs ${num_epochs} \
    --trials ${trials} --train_mode ${train_mode} --seed 123456

# check example output
cd ${example_dir}
python -u CBIG_DELSSOME_check_CMAES_pFIC_example_results.py

conda deactivate
