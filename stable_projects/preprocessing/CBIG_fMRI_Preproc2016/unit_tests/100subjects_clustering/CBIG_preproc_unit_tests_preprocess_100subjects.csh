#! /bin/csh -f
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set outdir = "/data/users/jingweil/storage/PreprocessingPipeline/release_test/20180104/100subjects_cluster/preproc_out"
set fmrinii_dir = "/mnt/eql/yeo3/data/GSP2016/CBIG_preproc_global_cen_bp/GSP_single_session/scripts/fmrinii"
set anat_dir = "/share/users/imganalysis/yeolab/data/GSP_release"
set config_file = "/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/100subjects_clustering/scripts/prepro.config"

set sub_list = "/mnt/eql/yeo3/data/GSP2016/CBIG_preproc_global_cen_bp/GSP_single_session/clustering/scripts/GSP_80_low_motion+20_w_censor.txt"

set curr_dir = `pwd`
set username = `whoami`
set work_dir = /data/users/$username/cluster/ 

echo $curr_dir
echo $username
echo $work_dir

if (! -e $work_dir) then
        mkdir -p $work_dir
endif

cd $work_dir


foreach curr_sub ("`cat $sub_list`")
	echo "curr_sub = $curr_sub"
	
	set cmd = "CBIG_preproc_fMRI_preprocess.csh -s $curr_sub -output_d $outdir -anat_s ${curr_sub}_FS -anat_d ${anat_dir} -fmrinii ${fmrinii_dir}/$curr_sub.fmrinii -config ${config_file}"
	echo $cmd | qsub -V -q circ-spool -l walltime=2:00:00,mem=4GB
	sleep 3s
	
end
