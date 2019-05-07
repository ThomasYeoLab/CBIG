#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <inDir> -d <idList> -s <scriptDir> -p <spmDir> [-q <queue>]
    Merge the grey matter images and create a mask.
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

echo 'Step 8: Merge and create mask.'

mkdir -p ${inDir}
logDir=${inDir}/logs/step8_merg_create_mask
mkdir -p ${logDir}
progressFile=${logDir}/progress.txt
> ${progressFile}

for id in `cat ${idList}`; do
    merge_list="${merge_list} ${inDir}/mri/mwp1${id}.nii"
    merge_list_smooth="${merge_list_smooth} ${inDir}/mri/s10_mwp1${id}.nii"
done

if [ -z "${queue}" ]; then
    fslmerge -t ${inDir}/mri/gm_merg ${merge_list}
    fslmerge -t ${inDir}/mri/gm_merg_s10 ${merge_list_smooth}

    fslmaths ${inDir}/mri/gm_merg -nan ${inDir}/mri/gm_merg
    fslmaths ${inDir}/mri/gm_merg_s10 -nan ${inDir}/mri/gm_merg_s10
    fslmaths ${inDir}/mri/gm_merg -Tmean -thr 0.1 -bin ${inDir}/mri/mask_fsl_0_1 -odt char
    echo "Done" >> ${progressFile}
else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N 'step8_merge_create_mask'
#PBS -l walltime=1:00:0
#PBS -l mem=16gb   
#PBS -e ${logDir}/merge.err
#PBS -o ${logDir}/merge.out
    
    cd ${scriptDir}

    fslmerge -t ${inDir}/mri/gm_merg ${merge_list}
    fslmerge -t ${inDir}/mri/gm_merg_s10 ${merge_list_smooth}

    fslmaths ${inDir}/mri/gm_merg -nan ${inDir}/mri/gm_merg
    fslmaths ${inDir}/mri/gm_merg_s10 -nan ${inDir}/mri/gm_merg_s10
    fslmaths ${inDir}/mri/gm_merg -Tmean -thr 0.1 -bin ${inDir}/mri/mask_fsl_0_1 -odt char
    echo "Done" >> ${progressFile}
EOJ
fi

./CBIG_MMLDA_waitUntilFinished.sh ${progressFile} 

echo 'Step 8: Merge and create mask. -- Finished.'

