#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -g <GMToNonlinTmpList> -o <outDir>
	- GMToNonlinTmpList	Text file with each line being the path to a GMToNonlinTmp image; e.g., ~/outputs/VBM/nonlinRegToNonlinTmp/GMToNonlinTmpList.txt
	- outDir		Output directory; e.g., ~/outputs/VBM/concatAndGenMask/
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":g:o:" opt; do
	case "${opt}" in
		g) GMToNonlinTmpList=${OPTARG};;
		o) outDir=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${GMToNonlinTmpList}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

echo 'Step 7: Concatenate Registered GM and Generate a Binary GM Mask'

mkdir -p ${outDir}
> ${outDir}GMToNonlinTmp_mod_4d_concatOrder.txt

# Prepare list
list=""
for GM in `cat ${GMToNonlinTmpList}`; do
	tmpVar="${GM##*/}" # discard everything before /
	filename="${tmpVar%_GMToNonlinTmp_mod.nii.gz}" # discard the suffix
	echo ${filename} >> ${outDir}GMToNonlinTmp_mod_4d_concatOrder.txt
	list="${list} ${GM}"
done

# Merge
fslmerge -t ${outDir}GMToNonlinTmp_mod_4d ${list}

# Generate binary GM mask
thr=0.05
fslmaths ${outDir}GMToNonlinTmp_mod_4d -Tmean -thr ${thr} -bin ${outDir}GMToNonlinTmp_mod_mean_binThr${thr} -odt char

echo 'Step 7: Concatenate Registered GM and Generate a Binary GM Mask -- Finished'

