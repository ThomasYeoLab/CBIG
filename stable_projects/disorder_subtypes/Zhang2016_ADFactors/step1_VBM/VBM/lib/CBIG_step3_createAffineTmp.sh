#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -g <GMToStdTmpList> -o <outDir>
	- GMToStdTmpList	Text file with each line being the path to a GMToStdTmp image; e.g., ~/outputs/VBM/affineRegToStdTmp/GMToStdTmpList.txt
	- outDir		Output directory; e.g., ~/outputs/VBM/affineTmp/
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":g:o:" opt; do
	case "${opt}" in
		g) GMToStdTmpList=${OPTARG};;
		o) outDir=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${GMToStdTmpList}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

echo 'Step 3: Create Study-Specific Affine Template'

mkdir -p ${outDir}
# Prepare list
tmpList=""
for img in `cat ${GMToStdTmpList}`; do
	tmpList="${tmpList} ${img}"
done

fslmerge -t ${outDir}GMToStdTmp_4d ${tmpList} # concatenate
fslmaths ${outDir}GMToStdTmp_4d -Tmean ${outDir}GMToStdTmp_4d_mean # average
fslswapdim ${outDir}GMToStdTmp_4d_mean -x y z ${outDir}GMToStdTmp_4d_mean_flipped # flip
fslmaths ${outDir}GMToStdTmp_4d_mean -add ${outDir}GMToStdTmp_4d_mean_flipped -div 2 ${outDir}affineTmp # re-average

echo 'Step 3: Create Study-Specific Affine Template -- Finished'

