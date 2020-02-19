#!/bin/sh

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

outDir=$1

# reference direcotry
workspace=${CBIG_TESTDATA_DIR}/stable_projects/\
disorder_subtypes/Sun2019_ADJointFactors/step1_SPM_VBM
refDir=${workspace}/results/VBM_use_given_template

imgList=${workspace}/data/ADNI2_bl_rawpath_set2.txt
idList=${workspace}/data/ADNI2_bl_rid_set2.txt
idListReorient=${workspace}/data/ADNI2_bl_rid_set2_reorient.txt
scriptDir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step1_SPM_VBM/code
spmDir=${CBIG_SPM_DIR}
newTemp="${workspace}/results/VBM_create_new_template/mri/wTemplate_1.nii,1"
reorientSample=${scriptDir}/Translation_Identity.mat
queue=circ-spool

# step1 convert from raw image to nifti file
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh -i ${imgList} -d ${idList} -s 1 \
-t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${outDir} -q ${queue}

# step1a reorient the image, you may need to repeat this step with different reorient sample
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh -i ${imgList} -d ${idList} -s 1a \
-t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${outDir} -r ${reorientSample} -q ${queue}

# create new id list with reorientation
${scriptDir}/CBIG_MMLDA_append_suffix_list.sh ${idList} _reorient ${idListReorient}

# step2-2a segmentation
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh -i ${imgList} -d ${idListReorient} \
-s 6_2a -t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${outDir} -q ${queue}

# after segmentation QC, you may need to update the idListReorient

# step3-12
${scriptDir}/CBIG_MMLDA_runVBM_givenTemp.sh -i ${imgList} -d ${idListReorient} \
-s 7_12 -t ${newTemp} -c ${scriptDir} -p ${spmDir} -o ${outDir} -q ${queue}

# compare results between your output and reference output
fslmaths ${outDir}/mri/gm_merg_MNI2mm_s10.nii.gz -sub ${refDir}/mri/gm_merg_MNI2mm_s10.nii.gz ${outDir}/diff.nii.gz
diff_out_ref=`fslstats ${outDir}/diff.nii.gz -m`
thre=0.00001
if (( $(echo "$diff_out_ref < $thre" | bc -l) )); then
    echo "[PASSED]"
else
    echo "[FAILDED]"
    echo "Mean difference between output and reference volumn: ${diff_out_ref}"
fi
