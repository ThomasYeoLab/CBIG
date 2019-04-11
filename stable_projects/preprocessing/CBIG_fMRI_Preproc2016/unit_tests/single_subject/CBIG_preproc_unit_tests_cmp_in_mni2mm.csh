#!/bin/csh
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set pipe_dir1 = "/mnt/eql/yeo1/CBIG_private_data/unit_tests/stable_projects/preprocessing"
set pipe_dir1 = "$pipe_dir1/CBIG_fMRI_Preproc2016/single_subject/data"
set pipe_name1 = "gt"
set pipe_stem1 = "_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_MNI2mm_sm6_finalmask"

set pipe_dir2 = $1
set pipe_name2 = "user-test"
set pipe_stem2 = $2

set subject_id = "Sub1116_Ses1"
set runs = "002 003"

set output_dir = $3

foreach run ($runs)

	set cmd = ( matlab -nodesktop -nosplash -r '"' 'addpath(fullfile(getenv('"'"CBIG_CODE_DIR"'"'), '"'"stable_projects"'"', '"'"preprocessing"'"', '"'"CBIG_fMRI_Preproc2016"'"', '"'"utilities"'"'))'; CBIG_preproc_compare_two_pipelines $pipe_dir1 $pipe_name1 $pipe_stem1 $pipe_dir2 $pipe_name2 $pipe_stem2 $subject_id $run $output_dir vol; exit; '"' )
	echo $cmd
	eval $cmd

end
