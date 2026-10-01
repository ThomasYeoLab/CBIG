#!/bin/bash
# This script submit a job to HPC for CBIG_fMRI_Preproc2016 motion correction unit test.
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

outdir=$1   # Your output directory
fmrinii_dir="$CBIG_TESTDATA_DIR/stable_projects/preprocessing"
fmrinii_dir="$fmrinii_dir/CBIG_fMRI_Preproc2016/single_subject/scripts/fmrinii"
config_file="$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016\
/unit_tests/single_subject/prepro_motion_correction.config"
anat_dir="$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data"

curr_sub="TMR00032_TMSD"

curr_dir=`pwd`
work_dir=$HOME/cluster/ 

echo $curr_dir
echo $work_dir

if [ ! -e $work_dir ] 
then
    mkdir -p $work_dir
fi

cd $work_dir

log_file="${outdir}/CBIG_preproc_unit_tests_call_fMRI_preproc_motion_correction.log"

cmd="source CBIG_init_conda; conda activate CBIG_py3;"
cmd="$cmd csh ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh"
cmd="$cmd -s $curr_sub -output_d $outdir -anat_s ${curr_sub}_FS -anat_d"
cmd="$cmd ${anat_dir} -fmrinii ${fmrinii_dir}/${curr_sub}.fmrinii -config ${config_file}"
cmd="$cmd -nocleanup | tee -a ${log_file} ;conda deactivate"

temp_script_file="${outdir}/temp_script.sh"
echo '#!/bin/bash' >> ${temp_script_file}
echo ${cmd} >> ${temp_script_file}
chmod 755 ${temp_script_file}

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script_file}" -walltime 5:00:00 -mem 32G \
-name "CBIG_preproc_unit_tests_mc"
