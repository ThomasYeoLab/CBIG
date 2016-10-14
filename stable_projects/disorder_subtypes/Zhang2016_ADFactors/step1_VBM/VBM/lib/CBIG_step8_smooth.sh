#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -v <4DVolume> -s <sigma> -o <outDir>
	- 4DVolume	Images to smooth (concatenated into a 4D volume); e.g., ~/outputs/VBM/concatAndGenMask/GMToNonlinTmp_mod_4d.nii.gz
	- sigma		Sigma (in mm) of the gaussian kernel for smoothing
	- outDir	Output directory; e.g., ~/outputs/VBM/smoothing/
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":v:s:o:" opt; do
	case "${opt}" in
		v) vol=${OPTARG};;
		s) sigma=${OPTARG};;
		o) outDir=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${vol}" ] || [ -z "${sigma}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

echo 'Step 8: Smooth GM Images'

mkdir ${outDir}

fslmaths ${vol} -s ${sigma} ${outDir}GMToNonlinTmp_mod_4d_s${sigma}

echo 'Step 8: Smooth GM Images -- Finished'

