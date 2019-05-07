#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <id_list> -s <script_dir> -o <out_dir> [-q <queue>]
    Register PET to T1 image.
    - id_list       Text file with each line being the id of a subject
    - script_dir    Directory of current script
    - out_dir       Output directory 
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:s:o:q:" opt; do
    case "${opt}" in
            i) id_list=${OPTARG};;
            s) script_dir=${OPTARG};;
            o) out_dir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${id_list}" ] || [ -z "${script_dir}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 4: Normalize PET image wrt cerebellum gm.'

step_name="step4_normalize_cerebellum_gm"

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

    PET_file=${sub_dir}/PET_in_T1.nii.gz
    mask_file=${sub_dir}/mask_T1_cerebellum_gm.nii.gz
    output_file=${sub_dir}/PET_in_T1_NormCere.nii.gz

    if [ -z "${queue}" ]; then
        export SUBJECTS_DIR=${out_dir}
        cd ${script_dir}

        # create cerebellum GM mask in native T1 space
        mri_binarize --i ${sub_dir}_FS/mri/aparc+aseg.mgz --match 8 47 \
        --o ${sub_dir}/mask_T1_cerebellum_gm.nii.gz >> ${log_file}
        # apply the mask to PET_in_T1 image and normalize the image by average cerebellum GM value
        matlab -nosplash -nodisplay -nodesktop -r \
        "clear;MRI_file='${PET_file}';Mask_file='${mask_file}'; \
        Output_file='${output_file}';CBIG_MMLDA_normalize_cerebellum_gm;exit;" \
        >> ${log_file}

        echo "${id}" >> ${progress_file}
    else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N "${step_name}"
#PBS -l walltime=1:00:0
#PBS -l mem=2gb   
#PBS -e ${sub_dir}/logs/${step_name}.err
#PBS -o ${sub_dir}/logs/${step_name}.out

    export SUBJECTS_DIR=${out_dir}
    cd ${script_dir}
    # create cerebellum GM mask in native T1 space
    mri_binarize --i ${sub_dir}_FS/mri/aparc+aseg.mgz --match 8 47 \
    --o ${sub_dir}/mask_T1_cerebellum_gm.nii.gz >> ${log_file}
    # apply the mask to PET_in_T1 image and normalize the image by average cerebellum GM value
    matlab -nosplash -nodisplay -nodesktop -r \
    "clear;MRI_file='${PET_file}';Mask_file='${mask_file}'; \
    Output_file='${output_file}';CBIG_MMLDA_normalize_cerebellum_gm;exit;" \
    >> ${log_file}

    echo "${id}" >> ${progress_file}
EOJ
    fi
done

# wait unit all jobs finished
total_num_job=`grep -c "^" ${id_list}`

${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/CBIG_MMLDA_wait_until_finished.sh \
${progress_file} ${total_num_job}

echo 'Step 4: Normalize PET image wrt cerebellum gm. -- Finished.'