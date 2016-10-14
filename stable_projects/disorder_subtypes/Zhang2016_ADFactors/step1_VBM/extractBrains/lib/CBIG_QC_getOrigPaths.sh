#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -l <brainNameList> -d <brainDir>
	- brainNameList		.txt file whose lines are brain filenames (e.g., 0002_bl_brain)
	- brainDir		Directory containing the BET outputs; e.g., ~/outputs/VBM/brains1/
Output .txt file is in the same folder as the input .txt file, with an additional suffix -- _orig.
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":l:d:" opt; do
	case "${opt}" in
		l) brainNameList=${OPTARG};;
		d) brainDir=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${brainNameList}" ] || [ -z "${brainDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

outFile=${brainNameList%".txt"}_orig.txt # assuming .txt
> ${outFile}
for brainName in `cat ${brainNameList}`; do
	origFilename=`cat ${brainDir}${brainName}_origPath.txt`
	echo ${origFilename} >> ${outFile}
done
