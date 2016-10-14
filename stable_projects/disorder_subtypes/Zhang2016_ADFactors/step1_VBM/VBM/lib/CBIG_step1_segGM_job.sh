#!/bin/sh

#PBS -l walltime=01:00:0
#PBS -l mem=6gb
# Will come into play only if this script is qsub'ed

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fast -R 0.3 -H 0.1 -o ${outDir}${filename} ${brain}

# Generate the partial GM list for template construction
if grep -Fxq "${brain}" ${brainList_tmp}; then
	echo ${outDir}${filename}_pve_1.nii.gz >> ${outDir}GMList_tmp.txt
fi

# Generate the full GM list for the next step (and indicate this job is done)
echo ${outDir}${filename}_pve_1.nii.gz >> ${outDir}GMList.txt
