#!/bin/csh
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set outdir = $1   # Your output directory
set fmrinii_dir = "$CBIG_TESTDATA_DIR/stable_projects/preprocessing"
set fmrinii_dir = "$fmrinii_dir/CBIG_fMRI_Preproc2016/single_subject/scripts/fmrinii"
set config_file = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016"
set config_file = "$config_file/unit_tests/single_subject/prepro_aCompCor.config"
set anat_dir = "$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data"

set curr_sub = "Sub1116_Ses1"

set curr_dir = `pwd`
set work_dir = $HOME/cluster/ 

echo $curr_dir
echo $work_dir

if (! -e $work_dir) then
    mkdir -p $work_dir
endif

cd $work_dir

set log_file = "${outdir}/CBIG_preproc_unit_tests_call_fMRI_preproc.log"

set cmd = "${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh"
set cmd = "$cmd -s $curr_sub -output_d $outdir -anat_s ${curr_sub}_FS -anat_d"
set cmd = "$cmd ${anat_dir} -fmrinii ${fmrinii_dir}/$curr_sub.fmrinii -config ${config_file} -nocleanup"
set cmd = "$cmd | tee -a ${log_file}"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 1:30:00 -mem 6G \
-name "CBIG_preproc_unit_tests_call_fMRI_preproc"
