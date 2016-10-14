#!/bin/sh

#PBS -l walltime=3:00:0
# Will come into play only if this script is qsub'ed

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fsl_reg ${GM} ${nonlinTmp} ${outDir}${filename}_GMToNonlinTmp -fnirt "--config=GM_2_MNI152GM_2mm.cnf --jout=${outDir}${filename}_JAC_nl"
# Each voxel of each registered GM image is multiplied by the Jacobian of the warp field
fslmaths ${outDir}${filename}_GMToNonlinTmp -mul ${outDir}${filename}_JAC_nl ${outDir}${filename}_GMToNonlinTmp_mod -odt float

# Generate GMToNonlinTmp_mod list for the next step (and indicate this job is done)
echo ${outDir}${filename}_GMToNonlinTmp_mod.nii.gz >> ${outDir}GMToNonlinTmpList.txt
