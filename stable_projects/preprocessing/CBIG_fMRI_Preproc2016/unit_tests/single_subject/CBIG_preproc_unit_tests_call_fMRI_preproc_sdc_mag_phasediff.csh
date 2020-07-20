#!/bin/csh
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set outdir = $1   # Your output directory
set fmrinii_dir = "$CBIG_TESTDATA_DIR/stable_projects/preprocessing"
set fmrinii_dir = "$fmrinii_dir/CBIG_fMRI_Preproc2016/single_subject/scripts/fmrinii"
set config_file = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016"
set config_file = "$config_file/unit_tests/single_subject/prepro_sdc_mag_phasediff.config"
set anat_dir = "$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data"

set curr_sub = "sub-NDARAA075AMK"

set curr_dir = `pwd`
set work_dir = $HOME/cluster/ 

echo $curr_dir
echo $work_dir

if (! -e $work_dir) then
    mkdir -p $work_dir
endif

cd $work_dir


set cmd = "CBIG_preproc_fMRI_preprocess.csh -s $curr_sub -output_d $outdir -anat_s ${curr_sub}_FS -anat_d"
set cmd = "$cmd ${anat_dir} -fmrinii ${fmrinii_dir}/${curr_sub}_task-rest.fmrinii -config ${config_file} -nocleanup"
echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l walltime=1:00:00,mem=6GB -m ae -N \
    CBIG_preproc_unit_tests_call_fMRI_preproc
