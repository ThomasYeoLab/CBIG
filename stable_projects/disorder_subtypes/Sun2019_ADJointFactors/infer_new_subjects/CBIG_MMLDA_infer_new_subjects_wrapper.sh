#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

./CBIG_MMLDA_infer_new_subjects.sh \
-i /data/users/nanbos/storage/CBIG_private/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/infer_new_subjects/input/T1_list.txt \
-d /data/users/nanbos/storage/CBIG_private/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/infer_new_subjects/input/id_list.txt \
-s /data/users/nanbos/storage/CBIG_private/stable_projects/disorder_subtypes/\
Sun2019_ADJointFactors/infer_new_subjects/input/subinfo.csv \
-k 3 \
-p /apps/arch/Linux_x86_64/spm/spm12 \
-o /data/users/nanbos/storage/test/infer_new \
-q circ-spool
