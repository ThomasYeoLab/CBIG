#!/bin/sh -f
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

sh ${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/bash_scripts/CBIG_resample_metric_to_164k.sh $1 $2
sh ${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/bash_scripts/CBIG_border_to_fsaverage.sh $1 $2
