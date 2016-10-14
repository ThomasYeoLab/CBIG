#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

brainList=/data/users/xzhang/storage/forPNASRelease/outputs/VBM_m24/brainList.txt
nuisanceVars=/data/users/xzhang/storage/forPNASRelease/outputs/VBM_m24/age_sex_icv_m24_matchBrainList.csv
K=3
outDir=/data/users/xzhang/storage/forPNASRelease/outputs_infNew/
queue=circ-spool
# Run
./CBIG_inferNew.sh -b ${brainList} -n ${nuisanceVars} -k ${K} -o ${outDir} -q ${queue}
