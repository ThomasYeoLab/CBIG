#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 <brainDir>
	- brainDir	Directory containing the BET outputs; e.g., ~/outputs/VBM/brains1/
Outputs are in \${brainDir}QC. Open \${brainDir}QC/index.html with a web browser to view a compiled version of all results.
" 1>&2; exit 1; }

# Reading in parameters
if [ -z "$1" ]; then
	echo Missing Parameters!
	usage
else
	brainDir=$1
fi

###########################################
# Main
###########################################

imgList=""
for brain in ${brainDir}*_brain.nii.gz; do
	# Get path of the original volume
	tmpVar="${brain##*/}" # discard everything before /
	name="${tmpVar%_brain.nii.gz}" # discard the extension
	orig=`cat ${brainDir}${name}_origPath.txt`
	# Concatenate to list
 	imgList="${imgList} ${orig} ${brain}"
done

slicesdir -o ${imgList}

mv ./slicesdir ${brainDir}QC
