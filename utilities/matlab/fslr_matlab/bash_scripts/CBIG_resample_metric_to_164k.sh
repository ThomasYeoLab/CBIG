#!/bin/bash
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set -x

ExtractDir=$CBIG_CODE_DIR/data/templates/surface/
SubjectDir=$1
SubjList=$2
ResamplingMethod=BARYCENTRIC

for Subject in $SubjList
do
    for Hemisphere in L R
    do
        MetricIn=$SubjectDir/"$Subject".$Hemisphere.32k_fs_LR.label.gii
        CurrentSphere=$ExtractDir/fs_LR_32k/cifti/standard.$Hemisphere.sphere.32k_fs_LR.surf.gii
        NewSphere=$ExtractDir/fs_LR_164k/cifti/standard.$Hemisphere.sphere.164k_fs_LR.surf.gii
        MetricOutput=$SubjectDir/"$Subject".$Hemisphere.164k_fs_LR.label.gii
        echo $MetricOutput
        wb_command -metric-resample \
        ${MetricIn} \
        $CurrentSphere \
        $NewSphere \
        $ResamplingMethod \
        $MetricOutput \
        -largest
    done
done
###################################################################
exit                
###################################################################
###################################################################
standard.L.sphere.164k_fs_LR.surf.gii
standard.L.sphere.32k_fs_LR.surf.gii
standard.R.sphere.164k_fs_LR.surf.gii
standard.R.sphere.32k_fs_LR.surf.gii

