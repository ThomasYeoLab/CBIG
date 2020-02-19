#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

out_dir=$1

workspace="${CBIG_TESTDATA_DIR}/stable_projects/\
disorder_subtypes/Sun2019_ADJointFactors/step3_PET_preprocess"
# reference directory
ref_dir=${workspace}/results

id_list=${workspace}/data/RID_Viscode_2sub.txt
T1_list=${workspace}/data/T1_path_2sub.txt
PET_list=${workspace}/data/PET_path_2sub.txt
deform_dir=${workspace}/data/
script_dir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step3_PET_preprocess
spm_dir=${CBIG_SPM_DIR}
queue=circ-spool

# run PET preprocessing
./CBIG_MMLDA_runPETpreprocess.sh -i ${id_list} -t ${T1_list} -p ${PET_list} \
-e 1_2 -d ${deform_dir} -s ${script_dir} -m ${spm_dir} -o ${out_dir} -q ${queue}

./CBIG_MMLDA_runPETpreprocess.sh -i ${id_list} -t ${T1_list} -p ${PET_list} \
-e 3_6 -d ${deform_dir} -s ${script_dir} -m ${spm_dir} -o ${out_dir} -q ${queue}

# compare results between your output and reference output
fslmaths ${out_dir}/merg_SPMVBMdeform.nii.gz -sub ${ref_dir}/merg_SPMVBMdeform.nii.gz ${out_dir}/diff.nii.gz
diff_out_ref=`fslstats ${out_dir}/diff.nii.gz -m`
thre=0.00001
if (( $(echo "$diff_out_ref < $thre" | bc -l) )); then
    echo "[PASSED]"
else
    echo "[FAILDED]"
    echo "Mean difference between output and reference volumn: ${diff_out_ref}"
fi
