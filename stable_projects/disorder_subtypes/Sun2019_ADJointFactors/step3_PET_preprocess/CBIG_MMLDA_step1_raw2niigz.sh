#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <id_list> -t <T1_list> -p <PET_list> -o <out_dir> [-q <queue>]
    Convert raw T1 and PET image to nifti file and rename it.
    - id_list       Text file with each line being the id of a subject
    - T1_list       Text file with each line being the path to a T1 image
    - PET_list      Text file with each line being the path to a PET image
    - out_dir       Output directory 
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:t:p:o:q:" opt; do
    case "${opt}" in
            i) id_list=${OPTARG};;
            t) T1_list=${OPTARG};;
            p) PET_list=${OPTARG};;
            o) out_dir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${id_list}" ] || [ -z "${T1_list}" ] || [ -z "${PET_list}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 1: Convert raw T1 image and PET image to nifti file.'

step_name="step1_raw2niigz"

# create empty progress file
progress_file=${out_dir}/progress/${step_name}.txt
mkdir -p ${out_dir}/progress
> ${progress_file}

# read the id from id list
i=0
for id in `cat ${id_list}`; do
    i=$(($i+1))
    id_array[$i]=${id}
done

# read the T1 scan list
i=0
for scan in `cat ${T1_list}`; do
    i=$(($i+1))
    T1_path[$i]=${scan}
done

# read the PET scan list
i=0
for scan in `cat ${PET_list}`; do
    i=$(($i+1))
    PET_path[$i]=${scan}
done

# convert raw image to nifti file
i=0
for id in `cat ${id_list}`; do
    i=$(($i+1))
    id=${id_array[$i]}
    sub_dir=${out_dir}/${id}
    log_file=${sub_dir}/logs/${step_name}.log
    mkdir -p ${sub_dir}/logs
    > ${log_file}

    if [ -z "${queue}" ]; then
        mri_convert ${T1_path[$i]} -odt float ${sub_dir}/T1.nii.gz >> ${log_file}
        mri_convert ${PET_path[$i]} -odt float ${sub_dir}/PET.nii.gz >> ${log_file}
        echo "${id}" >> ${progress_file}
    else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N "${step_name}"
#PBS -l walltime=1:00:0
#PBS -l mem=2gb   
#PBS -e ${sub_dir}/logs/${step_name}.err
#PBS -o ${sub_dir}/logs/${step_name}.out

    mri_convert ${T1_path[$i]} -odt float ${sub_dir}/T1.nii.gz >> ${log_file}
    mri_convert ${PET_path[$i]} -odt float ${sub_dir}/PET.nii.gz >> ${log_file}
    echo "${id}" >> ${progress_file}
EOJ
    fi
done

# wait unit all jobs finished
total_num_job=`grep -c "^" ${id_list}`

${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/CBIG_MMLDA_wait_until_finished.sh \
${progress_file} ${total_num_job}

echo 'Step 1: Convert raw T1 image and PET image to nifti file. -- Finished.'