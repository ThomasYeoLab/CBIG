#!/bin/sh

#PBS -l walltime=01:00:0
# Will come into play only if this script is qsub'ed

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

stdTmp=${FSL_DIR}/data/standard/tissuepriors/avg152T1_gray
fsl_reg ${GM} ${stdTmp} ${outDir}${filename}_GMToStdTmp -a # Register to template

# Generate GMToStdTemp list for the next step (and indicate this job is done)
echo ${outDir}${filename}_GMToStdTmp.nii.gz >> ${outDir}GMToStdTmpList.txt

