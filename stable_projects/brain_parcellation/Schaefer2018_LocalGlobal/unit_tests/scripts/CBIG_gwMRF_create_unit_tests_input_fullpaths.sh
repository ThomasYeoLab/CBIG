#!/bin/bash

# This script creates a list that contains the full paths to unit test subjects's surf data. 
# Each line is one subject with different runs.
# Written by Xiaoxuan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



# set paths
data_dir=$CBIG_TESTDATA_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/data
output_dir=${1}

# create a file to store full paths
#output_file="${output_dir}/unit_tests_input_fullpaths.csv"
output_file="${output_dir}/unit_tests_input_fullpaths_short.csv"

if [ -f ${output_file} ]; then
    rm ${output_file}
fi
touch ${output_file}



# loop through all subjects, write full paths of each run in output file
cd ${data_dir}
subject_list=$( ls | grep Sub )

#for subject in ${subject_list}

#(to match the order of alex's first 10 subjects)
#for subject in "Sub0734_Ses1" "Sub1116_Ses1" "Sub1361_Ses1" "Sub0538_Ses1" "Sub0654_Ses1" 
#"Sub0207_Ses1" "Sub1537_Ses1" "Sub0282_Ses1" "Sub0390_Ses1" "Sub0458_Ses1"
for subject in "Sub0734_Ses1" "Sub0207_Ses1"
do

    data_stem_002="lh.${subject}_bld002_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6.nii.gz"
    data_stem_003="lh.${subject}_bld003_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6.nii.gz"

    full_path_002="${data_dir}/${subject}/surf/${data_stem_002}"
    full_path_003="${data_dir}/${subject}/surf/${data_stem_003}"

    if [ -f "${full_path_002}" ] && [ -f "${full_path_003}" ]; then
        echo -n "${full_path_002} " >> ${output_file}
        echo -e "${full_path_003}" >> ${output_file}
    fi

    if [ -f "${full_path_002}" ] && [ ! -f "${full_path_003}" ]; then
        echo -e "${full_path_002}" >> ${output_file}
    fi

    if [ ! -f "${full_path_002}" ] && [ -f "${full_path_003}" ]; then
        echo -e "${full_path_003}" >> ${output_file}
    fi

done

