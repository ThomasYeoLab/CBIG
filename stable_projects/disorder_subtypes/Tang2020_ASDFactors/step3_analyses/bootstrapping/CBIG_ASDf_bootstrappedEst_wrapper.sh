#!/bin/bash
# 
# Wrapper script to run polarLDA estimate on bootstrapped samples
# 
# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Take input variables
input_dir=$1 # absolute directory where the docs are saved 
out_dir=$2
K=$3 # number of factors
N=$4 # number of resamples
cluster=$5 # cluster name

proj_dir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Tang2020_ASDFactors
code_dir=${proj_dir}/step2_polarLDA
model_name=${proj_dir}/data_release/files_polarLDA/final
inf_settings=${code_dir}/CBIG_ASDf_polarLDA_infSettings.txt

for (( i=1; i<=${N}; i++ ))
do
    corpusDir_est=${input_dir}/resampled_${i}/dx1.dat
    output_dir=${out_dir}/resampled_${i}
    mkdir -p ${output_dir}
    output_dir_step2a=${output_dir}/estimate

    sh ${code_dir}/CBIG_ASDf_polarLDA_est_initFromModel.sh \
    -d ${corpusDir_est} \
    -t ${inf_settings} \
    -k ${K} \
    -i model \
    -m ${model_name} \
    -p ${code_dir} \
    -o ${output_dir_step2a} \
    -q ${cluster}

done
