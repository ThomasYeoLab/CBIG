#!/bin/sh -f
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

sh ${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/bash_scripts/CBIG_border_to_fs_lr.sh $1 $2 $3 #project from freesurfer to fslr 164k
sh ${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/bash_scripts/CBIG_resample_metric.sh $1 $2 #resample to 32k fslr
