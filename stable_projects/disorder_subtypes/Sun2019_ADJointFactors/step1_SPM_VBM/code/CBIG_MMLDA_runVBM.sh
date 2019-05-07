#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <imgList> -d <idList> -s <step> -c <scriptDir> -p <spmDir> -o <outDir> [-r <reorient>] [-q <queue>]
    Run VBM with 4 stages. For the 1st stage, it will run step1 and convert raw image to nifti file.
    For the 2nd stage, it will run step1a and the user need to reorient the images until all images have 
    origins in AC. And the author need to update the idList by adding suffix '_reorient' to the original id using 
    CBIG_MMLDA_append_suffix_list.sh. For the 3rd stage, it will run step 2, 2a
    it runs initial segmentation and the user should check the segmentation results 
    by ${outDir}/slicesdir/index.html. After checking, the user should update the 
    idList by removing the id of failed subjects. For the 4th stage, it runs remaining 3 to 12 steps.

    - imgList           Text file with each line being the path to a T1 image 
    - idList            Text file with each line being the ID to a subject. The 
                        id list should correspond to the img list line by line. 
    - step              Run step 1, 1a, 1b, 2_2a, 3_12. If the user wants to run a single
                        step, just run the code of that step in the wrapper.
    - scriptDir         Directory of current script
    - spmDir            Director of spm software
    - outDir            Output directory; e.g., ~/outputs/VBM/
    - reorient          If step = 1a, reorientation sample (by following youtube video) used for reorientation
                        a mat file similar to Translation_Identity.mat
                        If step = 1b, reorientation matrix list, each line is a path of reorientation matrix.
    - queue             (Optional) if you have a cluster, use it to specify the 
                        queue to which you want to qsub these jobs; if not provided,
                        jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:d:s:c:p:o:r:q:" opt; do
    case "${opt}" in
            i) imgList=${OPTARG};;
            d) idList=${OPTARG};;
            s) step=${OPTARG};;
            c) scriptDir=${OPTARG};;
            p) spmDir=${OPTARG};;
            o) outDir=${OPTARG};;
            r) reorient=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${imgList}" ] || [ -z "${idList}" ] || [ -z "${step}" ] || \
    [ -z "${scriptDir}" ] || [ -z "${spmDir}" ] || [ -z "${outDir}" ]; then
    echo Missing Parameters!
    usage
fi
if [ -z "${queue}" ]; then
    queue=""
fi

###########################################
# Main
###########################################
cd ${scriptDir}

if [ ${step} == "1" ]; then
    ###### Step 1: Convert raw image to nifti file 
    ./CBIG_MMLDA_step1_raw2nii.sh -i ${imgList} -d ${idList} -o ${outDir} -q ${queue}

elif [ ${step} = "1a" ]; then
    ###### Step 1a: Reorient the T1 image 
    if [ ! -z ${reorient} ]; then
        ./CBIG_MMLDA_step1a_reorient.sh -i ${outDir} -d ${idList} -r ${reorient} \
        -s ${scriptDir} -p ${spmDir} -q ${queue}
    fi

elif [ ${step} = "1b" ]; then
    ###### Step 1b: Apply reorientation matrix to T1 image
    if [ ! -z ${reorient} ]; then
        ./CBIG_MMLDA_step1b_apply_reorient_matrix.sh -i ${outDir} -d ${idList} \
        -r ${reorient} -s ${scriptDir} -p ${spmDir} -q ${queue}
    fi

elif [ ${step} == "2_2a" ]; then
    ###### Step 2: First segmentation using default templates. 
    ./CBIG_MMLDA_step2_segGM.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 2a: Check segmentation results.
    ./CBIG_MMLDA_step2a_check_segment.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

elif [ ${step} == "3_12" ]; then
    ###### Step 3: Create customized DARTEL template. 
    ./CBIG_MMLDA_step3_dartel.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 4: Population to ICBM. 
    ./CBIG_MMLDA_step4_population2ICBM.sh -i ${outDir} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 5: Deformation. 
    ./CBIG_MMLDA_step5_deformation.sh -i ${outDir} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 6: Segmentation using new templates. 
    newTemp=${outDir}"/mri/wTemplate_1.nii,1"
    ./CBIG_MMLDA_step6_seg_new_template.sh -i ${outDir} -d ${idList} -t ${newTemp} \
    -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 7: Smoothing. 
    ./CBIG_MMLDA_step7_smooth.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 8: Merge and create mask. 
    ./CBIG_MMLDA_step8_merge_create_mask.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 9: Compute GM and ICV.
    ./CBIG_MMLDA_step9_compute_GM_ICV.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 10: Downsample to MNI 2mm.
    ./CBIG_MMLDA_step10_downsample_2mm.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 11: Smoothing 2mm. 
    ./CBIG_MMLDA_step11_smooth_2mm.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}

    ###### Step 12: Merge and create mask 2mm. 
    ./CBIG_MMLDA_step12_merge_create_mask_2mm.sh -i ${outDir} -d ${idList} -s ${scriptDir} -p ${spmDir} -q ${queue}
fi