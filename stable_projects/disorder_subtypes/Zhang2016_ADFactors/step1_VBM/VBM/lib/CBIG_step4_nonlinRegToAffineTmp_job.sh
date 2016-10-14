#!/bin/sh

#PBS -l walltime=3:00:0
# Will come into play only if this script is qsub'ed

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fsl_reg ${GM} ${affineTmp} ${outDir}${filename}_GMToAffineTmp -fnirt "--config=GM_2_MNI152GM_2mm.cnf"

# Generate GMToAffineTmp list for the next step (and indicate this job is done)
echo ${outDir}${filename}_GMToAffineTmp.nii.gz >> ${outDir}GMToAffineTmpList.txt

