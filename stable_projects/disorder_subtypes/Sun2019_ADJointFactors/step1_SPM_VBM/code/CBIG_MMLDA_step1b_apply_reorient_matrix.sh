#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <inDir> -d <idList> -r <reorientMatrixList> -s <scriptDir> -p <spmDir> [-q <queue>]
    Apply existing reorient matrix to to T1 images.
    - inDir                 Input directory of T1 images
    - idList                Text file with each line being the id of a subject
    - reorientMatrixList    Reorientation matrix list. Each line is a file path of a mat file.
                            It's similar to Translation_Identity.mat
    - scriptDir             Script directory of current script
    - spmDir                SPM installation directory
    - queue                 (Optional) if you have a cluster, use it to specify the 
                            queue to which you want to qsub these jobs; if not provided,
                            jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:d:r:s:p:q:" opt; do
    case "${opt}" in
            i) inDir=${OPTARG};;
            d) idList=${OPTARG};;
            r) reorientMatrixList=${OPTARG};;
            s) scriptDir=${OPTARG};;
            p) spmDir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${inDir}" ] || [ -z "${idList}" ] || [ -z "${reorientMatrixList}" ] || \
    [ -z "${scriptDir}" ] || [ -z "${spmDir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo 'Step 1b: Apply reorientation matrix.'

mkdir -p ${inDir}
logDir=${inDir}/logs/step1b_apply_reorient_matrix
mkdir -p ${logDir}
progressFile=${logDir}/progress.txt
> ${progressFile}

# read the id from id list
i=0
for id in `cat ${idList}`; do
    i=$(($i+1))
    id_array[$i]=${id}
done

# read the reorient matrix list
i=0
for mat in `cat ${reorientMatrixList}`; do
    i=$(($i+1))
    reorient_matrix_array[$i]=${mat}
done

i=0
for id in `cat ${idList}`; do
    i=$(($i+1))
    logFile=${logDir}/${id}.log
    imgPath=${inDir}/${id}.nii
    imgPathReorient=${inDir}/${id}_reorient.nii
    reorient_matrix=${reorient_matrix_array[$i]}
    cp ${imgPath} ${imgPathReorient}

    if [ -z "${queue}" ]; then
        matlab -nodisplay -nosplash -r \
        "spm_dir='${spmDir}';script_dir='${scriptDir}';\
        image_path='${imgPathReorient}';reorient_matrix='${reorient_matrix}';\
        CBIG_MMLDA_apply_reorient_matrix;exit;" > ${logFile}
        echo "${id}" >> ${progressFile}
    else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N 'step1b_apply_reorient_matrix'
#PBS -l walltime=1:00:0
#PBS -l mem=2gb   
#PBS -e ${logDir}/${id}.err
#PBS -o ${logDir}/${id}.out
    
    cd ${scriptDir}
    
    matlab -nodisplay -nosplash -r \
    "spm_dir='${spmDir}';script_dir='${scriptDir}';\
    image_path='${imgPathReorient}';reorient_matrix='${reorient_matrix}';\
    CBIG_MMLDA_apply_reorient_matrix;exit;" > ${logFile}
    echo "${id}" >> ${progressFile}
EOJ
    fi
done

./CBIG_MMLDA_waitUntilFinished.sh ${progressFile} ${idList}

echo 'Step 1b: Apply reorientation matrix. -- Finished.'