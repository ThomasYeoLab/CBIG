#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <inDir> -d <idList> -s <scriptDir> -p <spmDir> [-q <queue>]
    Downsample the grey matter from 1.5mm to 2mm.
    - inDir         Input directory of T1 images
    - idList        Text file with each line being the id of a subject
    - scriptDir     Script directory of current script
    - spmDir        SPM installation directory
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:d:s:p:q:" opt; do
    case "${opt}" in
            i) inDir=${OPTARG};;
            d) idList=${OPTARG};;
            s) scriptDir=${OPTARG};;
            p) spmDir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${inDir}" ] || [ -z "${idList}" ] || [ -z "${scriptDir}" ] || [ -z "${spmDir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 10: Downsample to MNI 2mm.'

mkdir -p ${inDir}
logDir=${inDir}/logs/step10_downsample
mkdir -p ${logDir}
progressFile=${logDir}/progress.txt
> ${progressFile}

for id in `cat ${idList}`; do
    logFile=${logDir}/${id}.log
    if [ -z "${queue}" ]; then
        ./CBIG_MMLDA_downsampleToMNI2mm.sh ${inDir}/mri/mwp1${id}.nii > ${logFile}
        echo "${id}" >> ${progressFile}
    else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N 'step10_downsample'
#PBS -l walltime=1:00:0
#PBS -l mem=8gb   
#PBS -e ${logDir}/${id}.err
#PBS -o ${logDir}/${id}.out
    
    cd ${scriptDir}
    
    ./CBIG_MMLDA_downsampleToMNI2mm.sh ${inDir}/mri/mwp1${id}.nii > ${logFile}
    echo "${id}" >> ${progressFile}
EOJ
    fi
done

./CBIG_MMLDA_waitUntilFinished.sh ${progressFile} ${idList}

echo 'Step 10: Downsample to MNI 2mm. -- Finished.'

