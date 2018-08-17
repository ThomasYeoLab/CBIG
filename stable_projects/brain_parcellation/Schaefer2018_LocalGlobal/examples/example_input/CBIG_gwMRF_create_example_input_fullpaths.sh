#!/bin/bash

# This simple script creates a list that contains the full paths to example subjects' surf data. Each line represents one subject with different runs.
# Assume the folder structure of CBIG respository is preserved
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



# Set paths
data_dir=${CBIG_CODE_DIR}/data/example_data/CoRR_HNU

output_dir=${1}
mkdir -p ${output_dir}


# Create a file to store example full paths
output_file="${output_dir}/example_input_fullpaths.csv"
touch ${output_file}


# Write full paths to example subjects' surf data, each line represents one subject with different runs
echo -n "${data_dir}/subj01/subj01_sess1/surf/lh.subj01_sess1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz " > ${output_file}
echo -e "${data_dir}/subj01/subj01_sess2/surf/lh.subj01_sess2_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz" >> ${output_file}

echo -n "${data_dir}/subj02/subj02_sess1/surf/lh.subj02_sess1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz " >> ${output_file}
echo -e "${data_dir}/subj02/subj02_sess2/surf/lh.subj02_sess2_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz" >> ${output_file}


