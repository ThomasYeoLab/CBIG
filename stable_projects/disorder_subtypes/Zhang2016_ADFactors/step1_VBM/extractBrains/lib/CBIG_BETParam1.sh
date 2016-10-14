#!/bin/sh

#PBS -l walltime=01:00:0
# Useful only when this script is qsub'ed

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Example BET procedure/parameters (to be optimized based on your data)
standard_space_roi ${img} ${outDir}${name}_roi.nii.gz -b
bet ${outDir}${name}_roi.nii.gz ${outDir}${name}_brain.nii.gz -f 0.3
