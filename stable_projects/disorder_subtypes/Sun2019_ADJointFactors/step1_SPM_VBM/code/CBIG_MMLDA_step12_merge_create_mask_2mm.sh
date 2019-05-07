#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <inDir> -d <idList> -s <scriptDir> -p <spmDir> [-q <queue>]
    Merge downsampled grey matter images and create a mask.
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

echo 'Step 12: Merge and create mask 2mm.'

mkdir -p ${inDir}
logDir=${inDir}/logs/step12_merg_create_mask_2mm
mkdir -p ${logDir}
progressFile=${logDir}/progress.txt
> ${progressFile}

merge_list=""
merge_list_smooth=""
for id in `cat ${idList}`; do
    merge_list="${merge_list} ${inDir}/mri/mwp1${id}_MNI2mm.nii"
    merge_list_smooth="${merge_list_smooth} ${inDir}/mri/s10_mwp1${id}_MNI2mm.nii"
done

if [ -z "${queue}" ]; then
    # merge 
    fslmerge -t ${inDir}/mri/gm_merg_MNI2mm ${merge_list}
    fslmerge -t ${inDir}/mri/gm_merg_MNI2mm_s10 ${merge_list_smooth}

    # create mask
    fslmaths ${inDir}/mri/gm_merg_MNI2mm -nan ${inDir}/mri/gm_merg_MNI2mm
    fslmaths ${inDir}/mri/gm_merg_MNI2mm_s10 -nan ${inDir}/mri/gm_merg_MNI2mm_s10
    fslmaths ${inDir}/mri/gm_merg_MNI2mm -Tmean -thr 0.1 -bin ${inDir}/mri/mask_fsl_0_1_MNI2mm -odt char

    # create brain template
    fslmaths ${spmDir}/toolbox/cat12/templates_1.50mm/Template_T1_IXI555_MNI152.nii \
             -mul ${spmDir}/toolbox/cat12/templates_1.50mm/brainmask.nii \
             ${inDir}/Template_T1_IXI555_MNI152_brain.nii
    ./CBIG_MMLDA_downsampleToMNI2mm.sh ${inDir}/Template_T1_IXI555_MNI152_brain.nii.gz

    # mask out the nonbrain part
    fslmaths ${inDir}/Template_T1_IXI555_MNI152_brain_MNI2mm.nii -thr 0.2 -bin ${inDir}/brainmask_MNI2mm_bin -odt char
    fslmaths ${inDir}/mri/mask_fsl_0_1_MNI2mm -mas ${inDir}/brainmask_MNI2mm_bin \
    ${inDir}/mri/mask_fsl_0_1_MNI2mm_maskbrain
    echo "Done" >> ${progressFile}
else
        qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N 'step12_merge_create_mask_2mm'
#PBS -l walltime=1:00:0
#PBS -l mem=16gb   
#PBS -e ${logDir}/merge.err
#PBS -o ${logDir}/merge.out
    
    cd ${scriptDir}

    # merge 
    fslmerge -t ${inDir}/mri/gm_merg_MNI2mm ${merge_list}
    fslmerge -t ${inDir}/mri/gm_merg_MNI2mm_s10 ${merge_list_smooth}

    # create mask
    fslmaths ${inDir}/mri/gm_merg_MNI2mm -nan ${inDir}/mri/gm_merg_MNI2mm
    fslmaths ${inDir}/mri/gm_merg_MNI2mm_s10 -nan ${inDir}/mri/gm_merg_MNI2mm_s10
    fslmaths ${inDir}/mri/gm_merg_MNI2mm -Tmean -thr 0.1 -bin ${inDir}/mri/mask_fsl_0_1_MNI2mm -odt char

    # create brain template
    fslmaths ${spmDir}/toolbox/cat12/templates_1.50mm/Template_T1_IXI555_MNI152.nii \
             -mul ${spmDir}/toolbox/cat12/templates_1.50mm/brainmask.nii \
             ${inDir}/Template_T1_IXI555_MNI152_brain.nii
    ./CBIG_MMLDA_downsampleToMNI2mm.sh ${inDir}/Template_T1_IXI555_MNI152_brain.nii.gz

    # mask out the nonbrain part
    fslmaths ${inDir}/Template_T1_IXI555_MNI152_brain_MNI2mm.nii -thr 0.2 -bin ${inDir}/brainmask_MNI2mm_bin -odt char
    fslmaths ${inDir}/mri/mask_fsl_0_1_MNI2mm -mas ${inDir}/brainmask_MNI2mm_bin \
    ${inDir}/mri/mask_fsl_0_1_MNI2mm_maskbrain
    echo "Done" >> ${progressFile}
EOJ
fi

./CBIG_MMLDA_waitUntilFinished.sh ${progressFile} 

echo 'Step 12: Merge and create mask 2mm. -- Finished.'

