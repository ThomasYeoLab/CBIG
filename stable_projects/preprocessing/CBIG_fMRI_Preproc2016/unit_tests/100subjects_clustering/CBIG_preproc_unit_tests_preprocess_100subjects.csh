#! /bin/csh -f
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set outdir = $1
set fmrinii_dir = ${outdir}/fmrinii
set anat_dir = "$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/100subjects_clustering"
set anat_dir = "$anat_dir/recon_all"
set config_file = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests"
set config_file = "${config_file}/100subjects_clustering/prepro.config"

set sub_list = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/100subjects_clustering"
set sub_list = "${sub_list}/GSP_80_low_motion+20_w_censor.txt"

set curr_dir = `pwd`
set work_dir = $HOME/cluster/ 

echo $curr_dir
echo $work_dir

if (! -e $work_dir) then
    mkdir -p $work_dir
endif

cd $work_dir

foreach curr_sub ("`cat $sub_list`")
	echo "curr_sub = $curr_sub"
	set log_file = "${outdir}/prep_100sub_ut_${curr_sub}.log"
	set cmd = "${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh"
	set cmd = "$cmd -s $curr_sub -output_d $outdir -anat_s ${curr_sub}_FS"
	set cmd = "$cmd -anat_d ${anat_dir} -fmrinii ${fmrinii_dir}/$curr_sub.fmrinii -config ${config_file}"
        set cmd = "$cmd | tee -a ${log_file}"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 3:00:00 -mem 4G -name "prep_100sub_ut" 
	sleep 3s
end
