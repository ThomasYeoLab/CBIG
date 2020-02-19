#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

./CBIG_MMLDA_infer_new_subjects.sh \
-i ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/infer_new_subjects/input/T1_list.txt \
-d ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/infer_new_subjects/input/id_list.txt \
-s ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/infer_new_subjects/input/subinfo.csv \
-k 3 \
-p ${CBIG_SPM_DIR} \
-o ~/storage/test/MMLDA/infer_new \
-q circ-spool
