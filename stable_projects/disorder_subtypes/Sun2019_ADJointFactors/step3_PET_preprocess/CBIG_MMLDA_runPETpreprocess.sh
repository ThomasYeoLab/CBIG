#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <id_list> -t <T1_list> -p <PET_list> -e <step> -d <deform_dir> \
    -s <script_dir> -m <spm_dir> -o <out_dir> [-q <queue>]
    Run PET preprocessing with 2 stages. The 1st stage will run step1 and step2,
    the user should check the reconall results. If some subjects fail the recon-all,
    the user need to exclude them from the id_list. The 2nd stage will run step3 to 
    step6.
    Step 1: Convert raw T1 image and PET image to nifti file
    Step 2: Freesurfer recon-all of T1 image
    Step 3: Freesurfer register PET to T1
    Step 4: Normalize PET image wrt cerebellum gm
    Step 5: Transform PET image in T1 space to MNI space
    Step 6: Merge PET image of subjects

    - id_list       Text file with each line being the id of a subject
    - T1_list       Text file with each line being the path to a T1 image
    - PET_list      Text file with each line being the path to a PET image
    - step          1_2 or 3_6
    - deform_dir    The directory of deformation field from SPM VBM
    - script_dir    The directory of current script
    - spm_dir       The directory of SPM12 software
    - out_dir       Output directory 
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:t:p:e:d:s:m:o:q:" opt; do
    case "${opt}" in
            i) id_list=${OPTARG};;
            t) T1_list=${OPTARG};;
            p) PET_list=${OPTARG};;
            e) step=${OPTARG};;
            d) deform_dir=${OPTARG};;
            s) script_dir=${OPTARG};;
            m) spm_dir=${OPTARG};;
            o) out_dir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${id_list}" ] || [ -z "${T1_list}" ] || [ -z "${PET_list}" ] || \
    [ -z "${step}" ] || [ -z "${deform_dir}" ] || [ -z "${script_dir}" ] || \
    [ -z "${spm_dir}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi
if [ -z "${queue}" ]; then
    queue=""
fi

###########################################
# Main
###########################################
cd ${script_dir}

if [ ${step} == "1_2" ]; then
    ###### Step 1: Convert raw T1 image and PET image to nifti file
    ./CBIG_MMLDA_step1_raw2niigz.sh -i ${id_list} -t ${T1_list} -p ${PET_list} -o ${out_dir} -q ${queue}

    ###### Step 2: Freesurfer recon-all of T1 image
    ./CBIG_MMLDA_step2_reconall.sh -i ${id_list} -o ${out_dir} -q ${queue}

elif [ ${step} == "3_6" ]; then
    ###### Step 3: Freesurfer register PET to T1
    ./CBIG_MMLDA_step3_PET2T1.sh -i ${id_list} -o ${out_dir} -q ${queue}

    ###### Step 4: Normalize PET image wrt cerebellum gm
    ./CBIG_MMLDA_step4_normalize_cerebellum_gm.sh -i ${id_list} -s ${script_dir} -o ${out_dir} -q ${queue}

    ###### Step 5: Transform PET image in T1 space to MNI space
    ./CBIG_MMLDA_step5_PET2MNI.sh -i ${id_list} -d ${deform_dir} -s ${script_dir} \
    -p ${spm_dir} -o ${out_dir} -q ${queue}

    ###### Step 6: Merge PET image of subjects
    ./CBIG_MMLDA_step6_merge.sh -i ${id_list} -o ${out_dir} -q ${queue}
fi
