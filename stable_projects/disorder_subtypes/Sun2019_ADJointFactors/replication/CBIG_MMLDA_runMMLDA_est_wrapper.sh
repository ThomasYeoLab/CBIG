#!/bin/bash
# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
doc_dir=$1
out_dir=$2

unit_test_path=${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors

if [ -z "$1" ]; then
    doc_dir=${unit_test_path}/step2_MMLDA/results/BrainBehavior2doc
fi
if [ -z "$2" ]; then
    out_dir=${unit_test_path}/step2_MMLDA/results/estimation
fi

###
# Run MMLDA estimation with different number of factors and random initializations
###
cd ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA

# ADNI2 bl AD k=2, 3, 4
for k in `seq 2 4`; do
    ./CBIG_MMLDA_est.sh \
    -a ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_brain.dat \
    -b ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_behavior.dat \
    -t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
    -k ${k} \
    -s 1 \
    -e 20 \
    -m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
    -o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1 \
    -q circ-spool
done

# ADNI2 bl AD k=3 10 fold
k=3
for i in `seq 1 10`; do
    ./CBIG_MMLDA_est.sh \
    -a ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_brain_train${i}.dat \
    -b ${doc_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_behavior_train${i}.dat \
    -t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
    -k ${k} \
    -s 1 \
    -e 20 \
    -m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
    -o ${out_dir}/ADNI2_bl_AD_meanCNstdALL_plus1_train${i} \
    -q circ-spool
done

# ADNI1 bl AD k=3
k=3
./CBIG_MMLDA_est.sh \
-a ${doc_dir}/ADNI1_bl_AD_meanCNstdALL_plus1_brain.dat \
-b ${doc_dir}/ADNI1_bl_AD_meanCNstdALL_plus1_behavior.dat \
-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt \
-k ${k} \
-s 1 \
-e 20 \
-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA \
-o ${out_dir}/ADNI1_bl_AD_meanCNstdALL_plus1 \
-q circ-spool
