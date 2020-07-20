#!/bin/csh
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set pipe_dir1 = "$CBIG_TESTDATA_DIR/stable_projects/preprocessing"
set pipe_dir1 = "$pipe_dir1/CBIG_fMRI_Preproc2016/single_subject/data"
set pipe_name1 = "gt"
set pipe_stem1 = "_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5"

set pipe_dir2 = $1
set pipe_name2 = "user-test"
set pipe_stem2 = $2

set subject_id = "sub-NDARBF851NH6"
set runs = "001 002"

set output_dir = $3

foreach run ($runs)

	set cmd = ( matlab -nodesktop -nosplash -r '"' 'addpath(fullfile(getenv('"'"CBIG_CODE_DIR"'"'), \
	    '"'"stable_projects"'"', '"'"preprocessing"'"', '"'"CBIG_fMRI_Preproc2016"'"', '"'"utilities"'"'))'; \
	    CBIG_preproc_compare_two_pipelines $pipe_dir1 $pipe_name1 $pipe_stem1 $pipe_dir2 $pipe_name2 $pipe_stem2 \
	    $subject_id $run $output_dir surf; exit; '"' )
	echo $cmd
	eval $cmd

end
