#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <id_list> -d <deform_dir> -s <script_dir> -p <spm_dir> -o <out_dir> [-q <queue>]
    Transform PET to MNI space by applying deformation field.
    - id_list       Text file with each line being the id of a subject
    - deform_dir    The directory of deformation field from SPM VBM
    - script_dir    The directory of current script
    - spm_dir       The directory of SPM12 software
    - out_dir       Output directory 
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:d:s:p:o:q:" opt; do
    case "${opt}" in
            i) id_list=${OPTARG};;
            d) deform_dir=${OPTARG};;
            s) script_dir=${OPTARG};;
            p) spm_dir=${OPTARG};;
            o) out_dir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${id_list}" ] || [ -z "${deform_dir}" ] || [ -z "${script_dir}" ] || \
    [ -z "${spm_dir}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 5: Transform PET image in T1 space to MNI space.'

step_name="step5_PET2MNI"

# create empty progress file
progress_file=${out_dir}/progress/${step_name}.txt
mkdir -p ${out_dir}/progress
> ${progress_file}

# convert raw image to nifti file
i=0
for id in `cat ${id_list}`; do

    sub_dir=${out_dir}/${id}
    log_file=${sub_dir}/logs/${step_name}.log
    mkdir -p ${sub_dir}/logs
    > ${log_file}

    PET_reorient_file=${sub_dir}/PET_in_T1_NormCere_reorientT1.nii
    deformation_field=${deform_dir}/y_${id}_reorient.nii

    if [ -z "${queue}" ]; then
        export SUBJECTS_DIR=${out_dir}
        cd ${script_dir}

        # reorient PET image 
        mri_vol2vol --mov ${sub_dir}/PET_in_T1_NormCere.nii.gz \
        --targ ${sub_dir}/T1.nii.gz \
        --o ${sub_dir}/PET_in_T1_NormCere_reorientT1.nii.gz \
        --regheader --no-save-reg >> ${log_file}
        # convert nii.gz to nii file
        mri_convert ${sub_dir}/PET_in_T1_NormCere_reorientT1.nii.gz \
        ${sub_dir}/PET_in_T1_NormCere_reorientT1.nii >> ${log_file}
        # transform PET image to MNI space by applying deformation field
        matlab -nosplash -nodisplay -nodesktop -r \
        "spm_dir='${spm_dir}';script_dir='${script_dir}';\
        image_path='${PET_reorient_file}';deformation_field='${deformation_field}';\
        CBIG_MMLDA_apply_deformation;exit;" >> ${log_file}
        # downsample from MNI1mm to MNI2mm
        mri_vol2vol --mov ${sub_dir}/wPET_in_T1_NormCere_reorientT1.nii \
        --s FSL_MNI152_FS4.5.0 \
        --targ ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz \
        --o ${sub_dir}/PET_in_T1_NormCere_SPMVBMdeform_MNI2mm.nii.gz \
        --regheader --no-save-reg >> ${log_file}

        echo "${id}" >> ${progress_file}
    else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N "${step_name}"
#PBS -l walltime=1:00:0
#PBS -l mem=8gb   
#PBS -e ${sub_dir}/logs/${step_name}.err
#PBS -o ${sub_dir}/logs/${step_name}.out

    export SUBJECTS_DIR=${out_dir}
    cd ${script_dir}

    # reorient PET image 
    mri_vol2vol --mov ${sub_dir}/PET_in_T1_NormCere.nii.gz \
    --targ ${sub_dir}/T1.nii.gz \
    --o ${sub_dir}/PET_in_T1_NormCere_reorientT1.nii.gz \
    --regheader --no-save-reg >> ${log_file}
    # convert nii.gz to nii file
    mri_convert ${sub_dir}/PET_in_T1_NormCere_reorientT1.nii.gz \
    ${sub_dir}/PET_in_T1_NormCere_reorientT1.nii >> ${log_file}
    # transform PET image to MNI space by applying deformation field
    matlab -nosplash -nodisplay -nodesktop -r \
    "spm_dir='${spm_dir}';script_dir='${script_dir}';\
    image_path='${PET_reorient_file}';deformation_field='${deformation_field}';\
    CBIG_MMLDA_apply_deformation;exit;" >> ${log_file}
    # downsample from MNI1mm to MNI2mm
    mri_vol2vol --mov ${sub_dir}/wPET_in_T1_NormCere_reorientT1.nii \
    --s FSL_MNI152_FS4.5.0 \
    --targ ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz \
    --o ${sub_dir}/PET_in_T1_NormCere_SPMVBMdeform_MNI2mm.nii.gz \
    --regheader --no-save-reg >> ${log_file}

    echo "${id}" >> ${progress_file}
EOJ
    fi
done

# wait unit all jobs finished
total_num_job=`grep -c "^" ${id_list}`

${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/CBIG_MMLDA_wait_until_finished.sh \
${progress_file} ${total_num_job}

echo 'Step 5: Transform PET image in T1 space to MNI space. -- Finished.'