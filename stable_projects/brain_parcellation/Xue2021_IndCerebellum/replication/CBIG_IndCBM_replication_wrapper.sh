#!/bin/sh
# CBIG_IndCBM_replication_wrapper.sh <output_dir>
# This function intends to replicate the published result in this study, which involves the individual cerebellar 
# parcellations for 2 subjects.
# Note that this eplication is only avaliable within CBIG lab.
# Input:
#   output_dir:  Path of output folder. 
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

output_dir=$1
mkdir -p ${output_dir}/logs
mkdir -p ${output_dir}/temp

# Get subject name from a list file
sub1=`sed -n '1p' $CBIG_IndCBM_REP_DIR/rep_sub_list.txt`
sub2=`sed -n '2p' $CBIG_IndCBM_REP_DIR/rep_sub_list.txt`
sub_list="${sub1} ${sub2}"
set_list="all discovery replication"

echo "==========Replication: Generating data lists=========="

for sub in ${sub_list}
do
    cmd="sh $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/replication/"
    cmd="${cmd}CBIG_IndCBM_generate_list.sh"
    cmd="${cmd} $CBIG_IndCBM_REP_DIR/Derivatives/${sub}/smoothing/smooth_4mm"
    cmd="${cmd} $CBIG_IndCBM_REP_DIR/Individuals/${sub}/PROC_fsaverage6"
    cmd="${cmd} ${output_dir}/List/${sub}"
    echo $cmd
    eval $cmd
done

for set in ${set_list}
do
    cmd="sh $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/replication/"
    cmd="${cmd}CBIG_IndCBM_create_MSHBM_list.sh"
    cmd="${cmd} ${output_dir}/List/${sub1}/${set} ${output_dir}/List/${sub2}/${set}"
    cmd="${cmd} ${output_dir}/MSHBM/${set}"
    echo $cmd
    eval $cmd
done

echo "==========Replication: Computing surface profile=========="

for set in ${set_list}
do
    job_log=${output_dir}/logs/surf_profile_${set}.log
    cmd="sh $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/"
    cmd="${cmd}CBIG_IndCBM_compute_profile.sh fsaverage6 ${output_dir}/MSHBM/${set} | tee -a ${job_log}"
    echo $cmd

    temp_script="${output_dir}/temp/temp_script_surf_profile_${set}.sh"
    if [ -f "${temp_script}" ]; then
        rm ${temp_script}
    fi
    echo '#!/bin/bash' >> ${temp_script}
    echo ${cmd} >> ${temp_script}
    chmod 755 ${temp_script}
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script}" -walltime 10:00:00 -mem 20G -name "IndCBM_profile"
done

# Wait till jobs are finished
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n 'IndCBM_profile'

echo "==========Replication: Generating MSHBM parameters=========="

lh_profile="${output_dir}/MSHBM/all/profiles/avg_profile/lh_fsaverage6_roifsaverage3_avg_profile.nii.gz"
rh_profile="${output_dir}/MSHBM/all/profiles/avg_profile/rh_fsaverage6_roifsaverage3_avg_profile.nii.gz"
if [ ! -d ${output_dir}/MSHBM/all/group ];then
    mkdir -p ${output_dir}/MSHBM/all/group
fi
matlab -nosplash -nodisplay -nodesktop -r " \
addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Xue2021_IndCerebellum')));\
load(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Xue2021_IndCerebellum', \
'replication', 'input', 'HCP_10networks.mat')); \
clustered=CBIG_IndCBM_generate_MSHBM_params('${lh_profile}', '${rh_profile}', lh_labels, rh_labels); \
save('${output_dir}/MSHBM/all/group/group.mat', 'lh_labels', 'rh_labels', 'clustered', 'colors'); \
exit;"

group_file=`realpath ${output_dir}/MSHBM/all/group/group.mat`
if [ ! -d ${output_dir}/MSHBM/discovery/group ];then
    mkdir -p ${output_dir}/MSHBM/discovery/group
fi
ln -s ${group_file} ${output_dir}/MSHBM/discovery/group/group.mat
if [ ! -d ${output_dir}/MSHBM/replication/group ];then
    mkdir -p ${output_dir}/MSHBM/replication/group
fi
ln -s ${group_file} ${output_dir}/MSHBM/replication/group/group.mat

echo "==========Replication: Generating surface parcellation=========="

matlab -nosplash -nodisplay -nodesktop -r " \
addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Xue2021_IndCerebellum')));\
addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', \
'Kong2019_MSHBM', 'step2_estimate_priors'))); \
addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', \
'Kong2019_MSHBM', 'lib'))); \
CBIG_MSHBM_estimate_group_priors('${output_dir}/MSHBM/all', 'fsaverage6', '2', '31', '10' , \
'conv_th', '1e-6', 'save_all', '1'); \
CBIG_MSHBM_estimate_group_priors('${output_dir}/MSHBM/discovery', 'fsaverage6', '2', '16', '10' , \
'conv_th', '1e-6', 'save_all', '1'); \
CBIG_MSHBM_estimate_group_priors('${output_dir}/MSHBM/replication', 'fsaverage6', '2', '15', '10' , \
'conv_th', '1e-6', 'save_all', '1'); \
CBIG_IndCBM_extract_MSHBM_result('${output_dir}/MSHBM/all'); \
CBIG_IndCBM_extract_MSHBM_result('${output_dir}/MSHBM/discovery'); \
CBIG_IndCBM_extract_MSHBM_result('${output_dir}/MSHBM/replication'); \
exit;"

echo "==========Replication: Create cerebellum templates=========="

cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum

for sub in ${sub_list}
do
    cmd="sh $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/"
    cmd="${cmd}CBIG_IndCBM_create_template.sh fsaverage6"
    cmd="${cmd} $CBIG_IndCBM_REP_DIR/Derivatives/${sub}/mask/${sub}_bin_mask.nii.gz ${output_dir}/template"
    echo $cmd
    eval $cmd
done

echo "==========Replication: Compute FC=========="

if [ ! -d "${output_dir}/FC" ];then
    mkdir -p ${output_dir}/FC
fi

for sub in ${sub_list}
do
    for set in ${set_list}
    do
        job_log=${output_dir}/logs/fc_${sub}_${set}.log

        lh_surf_file="${output_dir}/List/${sub}/${set}/lh_${set}_list.txt"
        rh_surf_file="${output_dir}/List/${sub}/${set}/rh_${set}_list.txt"
        vol_file="${output_dir}/List/${sub}/${set}/MNI_${set}_list.txt"
        template="${output_dir}/template/${sub}_bin_mask_fsaverage6_cerebellum_template.dscalar.nii"
        output_file="${output_dir}/FC/${sub}_${set}_vol2surf"

        cmd="matlab -nosplash -nodisplay -nodesktop -r \" "
        cmd="${cmd} addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', "
        cmd="${cmd} 'Xue2021_IndCerebellum'))); "
        cmd="${cmd} CBIG_IndCBM_compute_vol2surf_fc('${lh_surf_file}', '${rh_surf_file}', '${vol_file}', "
        cmd="${cmd} '${template}', '${output_file}'); exit; \" | tee -a ${job_log}"
        echo $cmd
        
        temp_script="${output_dir}/temp/temp_script_fc_${sub}_${set}.sh"
        if [ -f "${temp_script}" ]; then
            rm ${temp_script}
        fi
        echo '#!/bin/bash' >> ${temp_script}
        echo ${cmd} >> ${temp_script}
        chmod 755 ${temp_script}
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script}" -walltime 30:00:00 -mem 250G -name "IndCBM_rep_fc"
    done
done

# Wait till jobs are finished
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n 'IndCBM_rep_fc'

echo "==========Replication: Cerebellar Parcellation=========="

if [ ! -d "${output_dir}/parcellation" ];then
    mkdir -p ${output_dir}/parcellation
fi

i=0
for sub in ${sub_list}
do
    i=$((i+1))
    for set in ${set_list}
    do
        job_log=${output_dir}/logs/wta_${sub}_${set}.log

        surf_label="${output_dir}/MSHBM/${set}/ind_parcellation/Ind_parcellation_MSHBM_sub${i}"
        vol2surf="${output_dir}/FC/${sub}_${set}_vol2surf"
        template="${output_dir}/template/${sub}_bin_mask_fsaverage6_cerebellum_template.dscalar.nii"
        par_dir="${output_dir}/parcellation"

        cmd="matlab -nosplash -nodisplay -nodesktop -r \" "
        cmd="${cmd} addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', "
        cmd="${cmd} 'Xue2021_IndCerebellum'))); "
        cmd="${cmd} CBIG_IndCBM_cerebellum_parcellation('${surf_label}', '${vol2surf}', '${template}', '${par_dir}', "
        cmd="${cmd} 'topX', '400', 'output_name', '${sub}_${set}_IndCBM_parcellation');  exit; \" | tee -a ${job_log}"
        echo $cmd

        temp_script="${output_dir}/temp/temp_script_wta_${sub}_${set}.sh"
        if [ -f "${temp_script}" ]; then
            rm ${temp_script}
        fi
        echo '#!/bin/bash' >> ${temp_script}
        echo ${cmd} >> ${temp_script}
        chmod 755 ${temp_script}
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script}" -walltime 2:00:00 -mem 3G -name "IndCBM_rep_wta"
    done
done

# Wait till jobs are finished
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n 'IndCBM_rep_wta'
rm -r "${output_dir}/temp"

echo "==========Replication: Done=========="
