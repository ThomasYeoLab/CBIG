#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -g <GMToAffineTmpList> -o <outDir>
	- GMToAffineTmpList	Text file with each line being the path to a GMToAffineTmp image; e.g., ~/outputs/VBM/nonlinRegToAffineTmp/GMToAffineTmpList.txt
	- outDir		Output directory; e.g., ~/outputs/VBM/nonlinTmp/
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":g:o:" opt; do
	case "${opt}" in
		g) GMToAffineTmpList=${OPTARG};;
		o) outDir=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${GMToAffineTmpList}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

echo 'Step 5: Create Study-Specific Nonlinear Template'

mkdir -p ${outDir}

# Prepare list
tmpList=""
for GM in `cat ${GMToAffineTmpList}`; do
	tmpList="${tmpList} ${GM}"
done

fslmerge -t ${outDir}GMToAffineTmp_4d ${tmpList} # concatenate
fslmaths ${outDir}GMToAffineTmp_4d -Tmean ${outDir}GMToAffineTmp_4d_mean # average
fslswapdim ${outDir}GMToAffineTmp_4d_mean -x y z ${outDir}GMToAffineTmp_4d_mean_flipped # flip
fslmaths ${outDir}GMToAffineTmp_4d_mean -add ${outDir}GMToAffineTmp_4d_mean_flipped -div 2 ${outDir}nonlinTmp # re-average

echo 'Step 5: Create Study-Specific Nonlinear Template -- Finished'

