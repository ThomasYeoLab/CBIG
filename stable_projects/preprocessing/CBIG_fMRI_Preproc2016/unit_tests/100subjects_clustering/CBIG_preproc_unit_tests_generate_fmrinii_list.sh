#!/bin/sh
# Written by Jingwei Li, Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

output_dir=$1
fmrinii_dir=${output_dir}/fmrinii

if [ ! -d $fmrinii_dir ]
then
    mkdir -p $fmrinii_dir
else
    rm -r $fmrinii_dir
    mkdir $fmrinii_dir
fi

subject_list=${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests
subject_list=${subject_list}/100subjects_clustering/GSP_80_low_motion+20_w_censor.txt

for subject in `cat ${subject_list}`
do
    if [ -e ${fmrinii_dir}/${subject}.fmrinii ]
    then
        rm ${fmrinii_dir}/${subject}.fmrinii
    fi

    cmd=${CBIG_TESTDATA_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016
    cmd=${cmd}/100subjects_clustering/preproc_out/${subject}/bold
    cd ${cmd}

    for i in 00?
    do
        echo "$i ${CBIG_TESTDATA_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/100subjects_clustering/preproc_out/${subject}/bold/$i/${subject}_bld${i}_rest.nii.gz" >> ${fmrinii_dir}/${subject}.fmrinii
    done
done

exit 0
