#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <id_list> -o <out_dir> [-q <queue>]
    Freesurfer recon-all of T1 image.
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

echo 'Step 6: Merge PET image of subjects.'

step_name="step6_merge"

# create empty progress file
progress_file=${out_dir}/progress/${step_name}.txt
mkdir -p ${out_dir}/progress
> ${progress_file}
log_file=${out_dir}/logs/${step_name}.log
mkdir -p ${out_dir}/logs

# create merge list
merge_list=""
for id in `cat ${id_list}`; do
    sub_dir=${out_dir}/${id}
    merge_list="${merge_list} ${sub_dir}/PET_in_T1_NormCere_SPMVBMdeform_MNI2mm.nii.gz"
done

if [ -z "${queue}" ]; then
    fslmerge -t ${out_dir}/merg_SPMVBMdeform ${merge_list}
    fslmaths ${out_dir}/merg_SPMVBMdeform -nan ${out_dir}/merg_SPMVBMdeform
    echo "Done" >> ${progress_file}
else
    qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N "${step_name}"
#PBS -l walltime=1:00:0
#PBS -l mem=8gb   
#PBS -e ${out_dir}/logs/${step_name}.err
#PBS -o ${out_dir}/logs/${step_name}.out

    fslmerge -t ${out_dir}/merg_SPMVBMdeform ${merge_list}
    fslmaths ${out_dir}/merg_SPMVBMdeform -nan ${out_dir}/merg_SPMVBMdeform
    echo "Done" >> ${progress_file}
EOJ
fi

# wait unit all jobs finished
${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/CBIG_MMLDA_wait_until_finished.sh \
${progress_file} 1

echo 'Step 6: Merge PET image of subjects. -- Finished.'