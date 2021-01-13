#!/bin/sh
# CBIG_IndCBM_compute_profile.sh <surf_mesh> <proj_dir>
# This function calls CBIG_MSHBM_generate_profiles.m and CBIG_MSHBM_avg_profiles.m to compute profiles based on a MSHBM
# folder structure. See CBIG_CODE_DIR/brain_parcellation/Kong2019_MSHBM for details.
# Input:
#   surf_mesh:   mesh name for your surface. Can be fsaverage5, fsaverage6, fsaverage, fs_LR_32k, fs_LR_164k.
#   proj_dir:    Path of MSHBM folder. This folder should include <proj_dir>/data_list/fMRI_list/?h.sub?_sess?.txt for
#                each hemisphere, each subject, each session. 
# Output:
#                Generated profile will be saved under <proj_dir>/profiles.
# Example:  CBIG_IndCBM_compute_profile.sh fsaverge5 <proj_dir>
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

surf_mesh=$1
proj_dir=$2
list_dir="${proj_dir}/data_list/fMRI_list"
output_dir="${proj_dir}/profiles"

sub_num=`ls ${list_dir} | grep lh_sub._sess1.txt | wc -l`
echo "${sub_num} subjects under ${proj_dir}"
for ((i=1;i<=$((sub_num));i++))
do
    max_sess_num=0;
    sess_num=`ls ${list_dir} | grep lh_sub${i} | wc -l`
    if [ ${sess_num} -gt $max_sess_num ];then
        max_sess_num=${sess_num}
    fi
    for ((j=1;j<=$((sess_num));j++))
    do
        if [ -d ${output_dir}/sub${i}/sess${j} ];then
            exist_file=`ls ${output_dir}/sub${i}/sess${j} | grep sub${i}_sess${j} | wc -l`
        else
            exist_file=0
        fi
        if [ $((exist_file)) == 2 ];then
            echo "Profiles for subject ${i}, session ${j} exist. Skip computing."
        else
            echo "Computing profile for subject ${i}, session ${j}..."
            matlab -nosplash -nodisplay -nodesktop -r " \
            addpath(genpath(fullfile('$CBIG_CODE_DIR', 'external_packages', 'SD'))); \
            addpath(genpath(fullfile('$CBIG_CODE_DIR', 'utilities', 'matlab'))); \
            addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', \
            'Kong2019_MSHBM', 'step1_generate_profiles_and_ini_params'))); \
            CBIG_MSHBM_generate_profiles('fsaverage3','${surf_mesh}', \
            '${proj_dir}','${i}','${j}','0'); \
            exit;"
        fi
    done
done
echo "Max session number is ${max_sess_num}."

# Compute average profile
lh_avg_profile="${output_dir}/avg_profile/lh_${surf_mesh}_roifsaverage3_avg_profile.nii.gz"
rh_avg_profile="${output_dir}/avg_profile/rh_${surf_mesh}_roifsaverage3_avg_profile.nii.gz"
if [ -f ${lh_avg_profile} -a -f ${rh_avg_profile} ];then
    echo "Average profile exist. Skip computing average profile."
else
    echo "Computing average profile..."
    matlab -nosplash -nodisplay -nodesktop -r " \
    addpath(genpath(fullfile('$CBIG_CODE_DIR', 'external_packages', 'SD'))); \
    addpath(genpath(fullfile('$CBIG_CODE_DIR', 'utilities', 'matlab'))); \
    addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', \
    'Kong2019_MSHBM', 'step1_generate_profiles_and_ini_params'))); \
    CBIG_MSHBM_avg_profiles 'fsaverage3' '${surf_mesh}' '${proj_dir}' '${sub_num}' '${sess_num}'; \
    exit;"
fi

# Make profile list
echo "Making profile list..."

profile_dir=${output_dir}
output_dir="${proj_dir}/profile_list/training_set"
if [ ! -d ${output_dir} ];then
    mkdir -p ${output_dir}
fi

hemilist="lh rh"
for hemi in ${hemilist}
do
    for ((i=1;i<$((max_sess_num))+1;i++))
    do
        output_file="${output_dir}/${hemi}_sess${i}.txt"
        if [ -f ${output_file} ];then
            rm -f ${output_file}
        fi
        for ((j=0;j<$((sub_num));j++))
        do
            profile="${profile_dir}/sub$((j+1))/sess${i}"
            profile="${profile}/${hemi}.sub$((j+1))_sess${i}_${surf_mesh}_roifsaverage3.surf2surf_profile.nii.gz"
            profile=`realpath ${profile}`
            if [ ! -f ${profile} ]; then
                echo "NONE" >> ${output_file}
            else
                echo ${profile} >> ${output_file}
            fi
        done
    done
done

echo "Done computing profile."

