#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <imgList> -d <idList> -o <outDir> [-q <queue>]
    Convert raw T1 image to nifti file and rename it with subject id.
    - imgList       Text file with each line being the path to a T1 image
    - idList        Text file with each line being the id of a subject
    - outDir        Output directory 
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:d:o:q:" opt; do
    case "${opt}" in
            i) imgList=${OPTARG};;
            d) idList=${OPTARG};;
            o) outDir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${imgList}" ] || [ -z "${idList}" ] || [ -z "${outDir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 1: Convert raw T1 image to nifti file.'

mkdir -p ${outDir}
logDir=${outDir}/logs/step1_raw2nii
mkdir -p ${logDir}
progressFile=${logDir}/progress.txt
> ${progressFile}

# read the id from id list
i=0
for id in `cat ${idList}`; do
    i=$(($i+1))
    id_array[$i]=${id}
done

# convert raw image to nifti file
i=0
for img in `cat ${imgList}`; do
    i=$(($i+1))
    id=${id_array[$i]}
    logFile=${logDir}/${id}.log
    if [ -z "${queue}" ]; then
        mri_convert ${img} -odt float ${outDir}/${id}.nii > ${logFile}
        echo "${outDir}/${id}.nii" >> ${progressFile}
    else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N 'step1_raw2nii'
#PBS -l walltime=1:00:0
#PBS -l mem=2gb   
#PBS -e ${logDir}/${id}.err
#PBS -o ${logDir}/${id}.out

    mri_convert ${img} -odt float ${outDir}/${id}.nii > ${logFile} 
    echo "${outDir}/${id}.nii" >> ${progressFile}
EOJ
    fi
done

./CBIG_MMLDA_waitUntilFinished.sh ${progressFile} ${idList}

echo 'Step 1: Convert raw T1 image to nifti file. -- Finished.'