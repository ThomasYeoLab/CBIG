#!/bin/bash
# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
doc_dir=$1
visualize_dir=$2
out_dir=$3

unit_test_path=${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors

if [ -z "$1" ]; then
    doc_dir=${unit_test_path}/step2_MMLDA/results/BrainBehavior2doc
fi
if [ -z "$2" ]; then
    visualize_dir=${unit_test_path}/step2_MMLDA/results/visualizeFactors
fi
if [ -z "$3" ]; then
    out_dir=${unit_test_path}/step2_MMLDA/results/inference
fi

###
# Run MMLDA inference based on ADNI2 bl k=3
###
cd ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA

k=3
rname=`ls ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/`

# ADNI2 bl AD k=3
./CBIG_MMLDA_inf.sh \
-a ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_brain.dat \
-b ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-d ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/${rname}/final \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1 \
-n ADNI2_bl_AD_meanCNstdALL_plus1 \
-q circ-spool

# ADNI2 bl AD k=3 10fold
for i in `seq 1 10`; do
    k=3
    rname=`ls ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_train${i}/k${k}/`
    ./CBIG_MMLDA_inf.sh \
    -a ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_brain_test${i}.dat \
    -b ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_behavior_test${i}.dat \
    -t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
    -k ${k} \
    -d ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_train${i}/k${k}/${rname}/final \
    -m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
    -o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_train${i} \
    -n ADNI2_bl_AD_meanCNstdALL_plus1_test${i} \
    -q circ-spool
done

k=3
rname=`ls ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/`

# ADNI2 bl MCI k=3
./CBIG_MMLDA_inf.sh \
-a ${doc_dir}/ADNI2_bl_MCI_meanCNstdALL_plus1_brain.dat \
-b ${doc_dir}/ADNI2_bl_MCI_meanCNstdALL_plus1_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-d ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/${rname}/final \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1 \
-n ADNI2_bl_MCI_meanCNstdALL_plus1 \
-q circ-spool

# ADNI2 bl CN k=3
./CBIG_MMLDA_inf.sh \
-a ${doc_dir}/ADNI2_bl_CN_meanCNstdALL_plus1_brain.dat \
-b ${doc_dir}/ADNI2_bl_CN_meanCNstdALL_plus1_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-d ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/${rname}/final \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1 \
-n ADNI2_bl_CN_meanCNstdALL_plus1 \
-q circ-spool

# ADNI2 m12 ALL
./CBIG_MMLDA_inf.sh \
-a ${doc_dir}/ADNI2_m12_ALL_meanCNstdALL_plus1_brain.dat \
-b ${doc_dir}/ADNI2_m12_ALL_meanCNstdALL_plus1_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-d ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/${rname}/final \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1 \
-n ADNI2_m12_ALL_meanCNstdALL_plus1 \
-q circ-spool

# ADNI2 PETtau ALL
./CBIG_MMLDA_inf.sh \
-a ${doc_dir}/ADNI2_PETtau_ALL_meanCNstdALL_plus1_brain.dat \
-b ${doc_dir}/ADNI2_PETtau_ALL_meanCNstdALL_plus1_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-d ${visualize_dir}/ADNI2_bl_AD_meanCNstdALL_plus1/k${k}/${rname}/final \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1 \
-n ADNI2_PETtau_ALL_meanCNstdALL_plus1 \
-q circ-spool
