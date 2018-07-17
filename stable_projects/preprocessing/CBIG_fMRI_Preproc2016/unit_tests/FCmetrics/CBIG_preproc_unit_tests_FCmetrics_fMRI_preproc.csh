#!/bin/csh
# Written by Jingwei Li, Yan Rui Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


#Sub0017 has 2 runs with no censor frame. Sub0735 has 1 run with censor frame. Sub1155 has 1 run with no censor frame. Sub1488 has 2 runs with censor frame



	set outdir = $1   # Your output directory
	set fmrinii_dir = "/mnt/eql/yeo3/data/GSP2016/CBIG_preproc_global_cen_bp/GSP_single_session/scripts/fmrinii"
	set config_file = "/data/users/tnyr/storage/Pre_Proc_Unit_Test/Preproc_Unit_Test_Pipeline.txt"
	set anat_dir = "/share/users/imganalysis/yeolab/data/GSP_release"

	

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

	set sub_array = (Sub0017_Ses1 Sub0735_Ses1 Sub1155_Ses1 Sub1488_Ses1)

	foreach i ($sub_array) 
	
		set curr_sub = $i

		set cmd = "CBIG_preproc_fMRI_preprocess.csh -s $curr_sub -output_d $outdir -anat_s ${curr_sub}_FS -anat_d ${anat_dir} -fmrinii ${fmrinii_dir}/$curr_sub.fmrinii -config ${config_file} -nocleanup"
		echo $cmd | qsub -V -q circ-spool -l walltime=02:00:00,mem=2GB

	end

	

	
