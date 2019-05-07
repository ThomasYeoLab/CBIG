#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <id_list> -o <out_dir> [-q <queue>]
    Register PET to T1 image.
    - id_list       Text file with each line being the id of a subject
    - out_dir       Output directory 
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:o:q:" opt; do
    case "${opt}" in
            i) id_list=${OPTARG};;
            o) out_dir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${id_list}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 3: Freesurfer register PET to T1.'

step_name="step3_PET2T1"

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

    if [ -z "${queue}" ]; then
        export SUBJECTS_DIR=${out_dir}
        mri_coreg --s ${id}_FS --mov ${sub_dir}/PET.nii.gz \
        --reg ${sub_dir}/PET2T1.reg.lta >> ${log_file}
        mri_vol2vol --lta ${sub_dir}/PET2T1.reg.lta --mov ${sub_dir}/PET.nii.gz \
        --fstarg --o ${sub_dir}/PET_in_T1.nii.gz >> ${log_file}
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
    mri_coreg --s ${id}_FS --mov ${sub_dir}/PET.nii.gz \
    --reg ${sub_dir}/PET2T1.reg.lta >> ${log_file}
    mri_vol2vol --lta ${sub_dir}/PET2T1.reg.lta --mov ${sub_dir}/PET.nii.gz \
    --fstarg --o ${sub_dir}/PET_in_T1.nii.gz >> ${log_file}
    echo "${id}" >> ${progress_file}
EOJ
    fi
done

# wait unit all jobs finished
total_num_job=`grep -c "^" ${id_list}`

${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/CBIG_MMLDA_wait_until_finished.sh \
${progress_file} ${total_num_job}

echo 'Step 3: Freesurfer register PET to T1. -- Finished.'