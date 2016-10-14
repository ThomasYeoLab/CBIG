#!/bin/csh -f 

# Author: Jingwei Li, Date: 2016/06/18
# example:
#     CBIG_general_cluster_fcMRI_surf2surf_profiles.csh -sd /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast -sub_ls /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast/scripts/GSP_newtestlist.txt -surf_stem fsaverage5 -n 7 -cluster_out GSP_100_low_motion_clusters

# This function is the main function that calls CBIG_create_subject_surf_list.csh, CBIG_compute_fcMRI_surf2surf_profiles_subjectlist.csh, and CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh in sequence. It assumes that the preprocessed surface data are located in sub_dir/subject/surf/. surf_stem is used to distinguish the data that you want to input from other surface data. 

set VERSION = '$Id: CBIG_general_cluster_surf2surf_profiles.csh v 1.0 2016/06/18 $'

set sub_dir = ""
set sub_list = ""
set surf_stem = ""
set scrub_flag = 0
set target = fsaverage5
set roi = fsaverage3
set num_clusters = ""
set out_dir = ""
set cluster_out = ""
set num_tries = 1000
set preproc_opt = "new"

set PrintHelp = 0;
if( $#argv == 0 ) goto usage_exit;
set n = `echo $argv | grep -e -help | wc -l`
if( $n != 0 ) then
	set PrintHelp = 1;
	goto usage_exit;
endif
set n = `echo $argv | grep -e -version | wc -l`
if( $n != 0 ) then
	echo $VERSION
	exit 0;
endif

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print os.path.realpath('$0')"`
set root_dir = `dirname $root_dir`

#########################################
# create log file
#########################################
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs
if( $scrub_flag == 0 ) then
	set LF = ${out_dir}/logs/CBIG_general_cluster_fcMRI_surf2surf_profiles_ncluster${num_clusters}_noscrub.log
else
	set LF = ${out_dir}/logs/CBIG_general_cluster_fcMRI_surf2surf_profiles_ncluster${num_clusters}_scrub.log
endif
rm -f $LF
touch $LF
echo "[Clustering]: logfile = $LF"


#########################################
# main function
#########################################
set cmd = "${root_dir}/CBIG_create_subject_surf_list.csh -sd ${sub_dir} -sub_ls ${sub_list} -surf_stem ${surf_stem} -out_dir ${out_dir} -preproc_opt ${preproc_opt}"
if( $scrub_flag == 1 ) then
	set cmd = "$cmd -outlier_stem ${outlier_stem}"
endif
echo $cmd |& tee -a $LF
eval $cmd |& tee -a $LF

set cmd = "${root_dir}/CBIG_compute_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list} -surf_ls ${out_dir}/lists/lh.surf${surf_stem}.list"
if( $scrub_flag == 1 ) then
	set cmd = "$cmd -outlier_ls ${out_dir}/lists/outlier${outlier_stem}.list"
endif
echo $cmd |& tee -a $LF
eval $cmd |& tee -a $LF

set cmd = "${root_dir}/CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list} -n ${num_clusters} -out ${cluster_out} -tries ${num_tries}"
if( $scrub_flag == 1 ) then
	set cmd = "$cmd -scrub_flag 1"
endif
echo $cmd |& tee -a $LF
eval $cmd |& tee -a $LF

exit 0


##########################################
# parse arguments
##########################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		# subjects directory (required)
		case "-sd":
			if( $#argv == 0 ) goto arg1err;
			set sub_dir = $argv[1]; shift;
			breaksw
			
		# subject list (required)
		case "-sub_ls":
			if( $#argv == 0 ) goto arg1err;
			set sub_list = $argv[1]; shift;
			breaksw
			
		# stem to identify input surface data name (required)
		case "-surf_stem":
			if( $#argv == 0 ) goto arg1err;
			set surf_stem = $argv[1]; shift;
			breaksw
			
		# outliers stem to specify outlier file name
		case "-outlier_stem":
			if( $#argv == 0 ) goto arg1err;
			set outlier_stem = $argv[1]; shift;
			set scrub_flag = 1;
			breaksw
			
		# target resolution (default is fsaverage5)
		case "-target":
			if( $#argv == 0 ) goto arg1err;
			set target = $argv[1]; shift;
			breaksw
			
		# ROI resolution (default is fsaverage3)
		case "-roi":
			if( $#argv == 0 ) goto arg1err;
			set roi = $argv[1]; shift;
			breaksw
			
		# number of clusters (required)
		case "-n":
			if( $#argv == 0 ) goto arg1err;
			set num_clusters = $argv[1]; shift;
			breaksw
			
		# output directory
		case "-out_dir":
			if( $#argv == 0 ) goto arg1err;
			set out_dir = $argv[1]; shift;
			breaksw
		
		# output filename of clusters (without extension, required)
		case "-cluster_out":
			if( $#argv == 0 ) goto arg1err;
			set cluster_out = $argv[1]; shift;
			breaksw
			
		case "-tries":
			if( $#argv == 0 ) goto arg1err
			set num_tries = $argv[1]; shift;
			breaksw
			
		# Assumption for preprocessing pipeline
		case "-preproc_opt":
			if( $#argv == 0 ) goto arg1err;
			set preproc_opt = $argv[1]; shift;
			breaksw
		
		default:
			echo "ERROR: flag $flag unrecognized."
			echo $cmdline
			exit 1
			breaksw	
		
	endsw
end

goto parse_args_return;


#####################################
# check parameters
#####################################
check_params:

if( $#sub_dir == 0 ) then
	echo "ERROR: subjects directory not specified."
	exit 1
endif

if( $#sub_list == 0 ) then
	echo "ERROR: subject list not specified."
	exit 1
endif

if( $#surf_stem == 0 ) then
	echo "ERROR: surface data stem not specified."
	exit 1
endif

if( $#num_clusters == 0 ) then
	echo "ERROR: number of clusters not specified."
	exit 1
endif

if( $#out_dir == 0 ) then
	echo "ERROR: output directory not specified."
	exit 1
endif

if( $#cluster_out == 0 ) then
	echo "ERROR: clustering output not specified."
	exit 1
endif

if( $#num_tries == 0 ) then
	echo "ERROR: number of tries not specified."
	exit 1
endif

goto check_params_return;


#####################################
# Error message
#####################################
arg1err:
  echo "ERROR: flag $flag requires one argument."
  exit 1


#####################################
# Usage exit
#####################################
usage_exit:

	echo ""
	echo "USAGE: CBIG_general_cluster_fcMRI_surf2surf_profiles.csh"
	echo ""
	echo "  Required arguments"
	echo "    -sd            sub_dir      : fMRI subject directory"
	echo "    -sub_ls        sub_list     : subject list"
	echo "    -surf_stem     surf_stem    : a stem that can identify the surface data that you want to use"
	echo "    -n             num_clusters : number of clusters"
	echo "    -out_dir       out_dir      : output directory"
	echo "    -cluster_out   cluster_out  : clustering output filename (full path, without extension)"
	echo ""
	echo "  Optional arguments"
	echo "    -tries         num_tries    : number of different random initialization (default is 1000)"
	echo "    -outlier_stem  outlier_stem : a stem that can identify the motion outliers file (without extension), e.g. '_FDRMS0.2_DVARS50_motion_outliers'"
	echo "    -target        target       : the resolution of clustering (default is fsaverage5)"
	echo "    -roi           roi          : the resolution of ROIs (defaule is fsaverage3)"
	echo "    -preproc_opt   preproc_opt  : assumption of preprocessing approach, choose from 'old' and 'new', 'old' means procsurffast file structure, 'new' means CBIG_fMRI_preprocess file structure. Default is 'new'"
	echo ""
	
	if ( $PrintHelp == 0 ) exit 1
	echo $VERSION
	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

  This function is the main function that calls CBIG_create_subject_surf_list.csh, CBIG_compute_fcMRI_surf2surf_profiles_subjectlist.csh, and CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh in sequence. The whole program implements the clustering algorithm of Yeo et al. 2011.
  It assumes that the preprocessed surface data are located in '${sub_dir}/${subject}/surf/' and motion outliers file is located in '${sub_dir}/${subject}/qc/'. If not, the user could create symbol links to these positions.
  -surf_stem is used to distinguish the data that you want to input from other surface data. 
  -outlier_stem is similar to -surf_stem. It is used to identify the motion outliers file. For example, for surface data preprocessed by 'CBIG_fMRI_preprocess.csh', the file name of motion outliers is '*_bld*_FDRMS0.2_DVARS50_motion_outliers.txt', then outlier_stem can be '_FDRMS0.2_DVARS50_motion_outliers'. This flag is optional. If it is specified, high motion frames (outliers) will be ignored when computing functional connectivity profiles.
  
