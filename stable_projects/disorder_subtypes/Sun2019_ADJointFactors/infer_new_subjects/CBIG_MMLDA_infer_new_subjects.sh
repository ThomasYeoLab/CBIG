#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <imgList> -d <idList> -s <subinfo> -k <k> -p <spmDir> -o <outDir> [-q <queue>]
    - imgList   Text file with each line being the path to a T1 image; e.g., input/T1_list.txt
    - idList    Text file with each line being the id of the subject; e.g. input/id_list.txt
    - subinfo   csv file which includes subject information, such as id, age, sex, behavior scores
                For the id column, it mush match the idList text file.
                For the format of the behavior scores, please check the input/subinfo.csv
                You have to fill in all behavior scores except the 'ADAS: Recall instruction, ANART, CFT: vegetable'
                columns. Because we do not use these three columns in following processing, you can just fill in these
                three columns with NaN.
    - k         Number of factors; possible values: 2, 3 and 4, because only these three models are released
    - spmDir    SPM directory; e.g., 
    - outDir    Output directory; e.g., ~/outputs/
    - queue     (Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs; if 
                not provided, jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:d:s:k:p:o:q:" opt; do
    case "${opt}" in
            i) imgList=${OPTARG};;
            d) idList=${OPTARG};;
            s) subinfo=${OPTARG};;
            k) k=${OPTARG};;
            p) spmDir=${OPTARG};;
            o) outDir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${imgList}" ] || [ -z "${idList}" ] || [ -z "${k}" ] || [ -z "${subinfo}" ] || \
    [ -z "${spmDir}" ] || [ -z "${outDir}" ]; then
    echo Missing Parameters!
    usage
fi
if [ "${k}" -ne 2 ] && [ "${k}" -ne 3 ] && [ "${k}" -ne 4 ]; then
    echo k can only be 2, 3 or 4!
    usage
fi

###########################################
# Main
###########################################

###### Step 1: VBM given GM Template
vbm_out_dir=${outDir}/vbm

newTemp=${CBIG_CODE_DIR}"/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/wTemplate_1.nii,1"
# newTemp="/mnt/eql/yeo11/data/ADNI_preprocess/Sun2019_SPMVBM/preprocessing/output/ADNI2_bl/mri/wTemplate_1.nii,1"
scriptDir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step1_SPM_VBM/code 

# convert from raw image to nifti file
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh \
-i ${imgList} -d ${idList} -s 1 -t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${vbm_out_dir} -q ${queue}
# segmentation
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh \
-i ${imgList} -d ${idList} -s 6_2a -t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${vbm_out_dir} -q ${queue}
# remaining steps
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh \
-i ${imgList} -d ${idList} -s 7_12 -t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${vbm_out_dir} -q ${queue}

###### Step 2: Converting VBM images and behavior scores to documents
doc_out_dir=${outDir}/doc
mkdir -p ${doc_out_dir}
logfile=${doc_out_dir}/CBIG_MMLDA_vbm_behavior_to_doc.log

gm_merge_path=${vbm_out_dir}/mri/gm_merg_MNI2mm_s10.nii.gz
gmvol_icv_path=${vbm_out_dir}/gmVolAndICV/id_gmVol_icv.mat
matlab -nodisplay -nosplash -r \
"subinfo_path='${subinfo}';gm_merge_path='${gm_merge_path}';\
gmvol_icv_path='${gmvol_icv_path}';doc_out_dir='${doc_out_dir}'; \
CBIG_MMLDA_vbm_behavior_to_doc;exit;" > ${logfile}

###### Step 3: MMLDA infer factor loadings given factors
mmlda_out_dir=${outDir}/mmlda_inf

cd ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA

./CBIG_MMLDA_inf.sh \
-a ${doc_out_dir}/doc_brain.dat \
-b ${doc_out_dir}/doc_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-d ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/\
ADNIDataRelease/MMLDA_files/ADNIGO2/k${k}/final \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${mmlda_out_dir} \
-n mmlda \
-q ${queue}

./CBIG_MMLDA_wait_until_finished.sh ${mmlda_out_dir}/k${k}_mmlda_progress.txt 1

###### Step 4: Normalize gamma file into probabilities
prob_out_dir=${outDir}/prob
mkdir -p ${prob_out_dir}
logfile=${prob_out_dir}/CBIG_MMLDA_gamma2prob.log

cd ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/infer_new_subjects

gamma_file=${mmlda_out_dir}/k${k}_inf_mmlda-gamma.dat
prob_file=${prob_out_dir}/k${k}_inf_mmlda_prob.csv
matlab -nodisplay -nosplash -r \
"gamma_file='${gamma_file}';prob_file='${prob_file}';\
CBIG_MMLDA_gamma2prob;exit;" > ${logfile}

gamma_file=${mmlda_out_dir}/k${k}_inf1_mmlda-gamma.dat
prob_file=${prob_out_dir}/k${k}_inf1_mmlda_prob.csv
matlab -nodisplay -nosplash -r \
"gamma_file='${gamma_file}';prob_file='${prob_file}';\
CBIG_MMLDA_gamma2prob;exit;" > ${logfile}

gamma_file=${mmlda_out_dir}/k${k}_inf2_mmlda-gamma.dat
prob_file=${prob_out_dir}/k${k}_inf2_mmlda_prob.csv
matlab -nodisplay -nosplash -r \
"gamma_file='${gamma_file}';prob_file='${prob_file}';\
CBIG_MMLDA_gamma2prob;exit;" > ${logfile}


