#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../lib/
smoothingSigma=4.246
queue=circ-spool


# VBM from scratch on 810 ADNI-1 baseline scans
brainList=/data/users/xzhang/storage/forPNASRelease/outputs/VBM_bl/brainList.txt
outDir=~/storage/forPNASRelease/outputs/VBM_bl/
./CBIG_runVBM.sh -b ${brainList} -s ${smoothingSigma} -o ${outDir} -q ${queue}


# VBM on 560 m24 follow-up scans using the template created above
brainList=/data/users/xzhang/storage/forPNASRelease/outputs/VBM_m24/brainList.txt
tmp=${outDir}nonlinTmp/nonlinTmp.nii.gz
outDir=~/storage/forPNASRelease/outputs/VBM_m24/
./CBIG_runVBM_givenTmp.sh -b ${brainList} -t ${tmp} -s ${smoothingSigma} -o ${outDir} -q ${queue}
