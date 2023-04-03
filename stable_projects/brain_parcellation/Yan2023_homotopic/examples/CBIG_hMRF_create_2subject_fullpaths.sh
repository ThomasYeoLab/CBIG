#!/bin/bash

# This script will generate the subject fullpaths for the two fake subjects: "Sub0734_Ses1" "Sub0207_Ses1"
# It checks through all possible runs (002, 003) of each subject if the corresponding fMRI file exists.
# A fullpath is then generated and tabulated into a csv file.
# The censor file paths will also be tabulated.
# Input arguments:
# 2) rh_out_file: directory to the left hemisphere output fullpath file
# 1) lh_out_file: directory to the right hemisphere output fullpath file
# 3) censor_file: directory to the censor file fullpath
#
# Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

lh_out_file=${1}
rh_out_file=${2}
censor_file=${3}

test_data_dir=${CBIG_TESTDATA_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/data
script_dir=`dirname $(readlink -f $0)`
input_data_dir=${script_dir}/input

if [[ -e ${lh_out_file} ]]; then
    echo "Output file already exists under output directory. File will be removed."
    rm -f "${lh_out_file}"
fi
touch ${lh_out_file}

if [[ -e ${rh_out_file} ]]; then
    echo "Output file already exists under output directory. File will be removed."
    rm -f "${rh_out_file}"
fi
touch ${rh_out_file}

if [[ -e ${censor_file} ]]; then
    echo "Output file already exists under output directory. File will be removed."
    rm -f "${censor_file}"
fi
touch ${censor_file}

#################
# Generate paths
#################
declare -a run_arr=("002" "003")
common_file_suffix='rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6'
for sub in "Sub0734_Ses1" "Sub0207_Ses1"; do
    #echo "Now is subject No. $sub"
    lh_cur_line=""
    rh_cur_line=""
    censor_cur_line=""
    for run in "${run_arr[@]}"; do
        lh_file_dir="${test_data_dir}/${sub}/surf/lh.${sub}_bld${run}_${common_file_suffix}.nii.gz"
        rh_file_dir="${test_data_dir}/${sub}/surf/rh.${sub}_bld${run}_${common_file_suffix}.nii.gz"
        temp_censor_file_dir="${input_data_dir}/fake_qc_file.txt"
        
        # this fake censoring file is shared across all runs
        if [[ -f "${lh_file_dir}" && -f "${rh_file_dir}" ]]; then
            lh_cur_line="${lh_cur_line} ${lh_file_dir}"
            rh_cur_line="${rh_cur_line} ${rh_file_dir}"
            censor_cur_line="${censor_cur_line} ${temp_censor_file_dir}" 
        fi
    done

    # write lh data
    if [[ -z ${lh_cur_line} ]]; then
        echo "Current subject ${sub} has no lh data."
    else
        echo "${lh_cur_line}" | xargs >> "${lh_out_file}"
    fi

    # write rh data
    if [[ -z ${rh_cur_line} ]]; then
        echo "Current subject ${sub} has no rh data."
    else
        echo "${rh_cur_line}" | xargs >> "${rh_out_file}"
    fi

    # write censoring file paths
    if [[ -z ${censor_cur_line} ]]; then
        echo "Current subject ${sub} has no fmri data."
    else
        echo "${censor_cur_line}" | xargs >> "${censor_file}"
    fi
done
