#!/bin/sh

#PBS -l walltime=01:00:0
# Useful only when this script is qsub'ed

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Example BET procedure/parameters (to be optimized based on your data)
bet ${img} ${outDir}${name}_brain.nii.gz -f 0.4 -B
