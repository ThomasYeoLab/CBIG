#!/bin/csh -f

# Author: Jingwei Li, Date: 2016/06/18

set VERSION = '$ID: CBIG_create_subject_surf_list.csh, v 1.0 2016/06/18'

set sub_dir = "";
set sub_list = "";
set fsaverage = "";
set surf_stem = ""
set outlier_stem = ""
set scrub_flag = 0;
set out_dir = ""
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


#################################
# output list configuration
#################################
set scripts_dir = "${out_dir}/lists"
mkdir -p $scripts_dir

set cmd = "CBIG_fMRI_create_data_list.csh -sd ${sub_dir} -sub_ls ${sub_list} -folder surf -data_stem ${surf_stem}.nii.gz -out_dir ${scripts_dir} -out_stem surf${surf_stem} -preproc_opt ${preproc_opt}"
echo $cmd
eval $cmd

if( $scrub_flag == 1 ) then
	set cmd = "CBIG_fMRI_create_data_list.csh -sd ${sub_dir} -sub_ls ${sub_list} -folder qc -data_stem ${outlier_stem}.txt -out_dir ${scripts_dir} -out_stem outlier${outlier_stem} -preproc_opt ${preproc_opt}"
	echo $cmd
	eval $cmd
endif



exit 0;


#################################
# parse arguments
#################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;

	switch($flag)
		# subject directory
		case "-sd":
			if( $#argv == 0 ) goto arg1err;
			set sub_dir = $argv[1]; shift;
			breaksw

		# subject list file
		case "-sub_ls":
			if( $#argv == 0 ) goto arg1err;
			set sub_list = $argv[1]; shift;
			breaksw

		# stem to identify input surface data name
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
			
		case "-out_dir":
			if( $#argv == 0 ) goto arg1err;
			set out_dir = $argv[1]; shift;
			breaksw
			
		# Assumption for preprocessing pipeline
		case "-preproc_opt":
			if( $#argv == 0 ) goto arg1err;
			set preproc_opt = $argv[1]; shift;
			breaksw

		default:
			echo "ERROR: Flag $flag unrecognized."
			echo $cmdline
			exit 1
			breaksw

	endsw
	
end

goto parse_args_return;


################################
# check passed parameters
################################
check_params:

if( $#sub_dir == 0 ) then
	echo "ERROR: subject directory not specified."
	exit 1;
endif

if( $#sub_list == 0 ) then
	echo "ERROR: subject list not specified."
	exit 1;
endif

if( $#surf_stem == 0 ) then
	echo "ERROR: surface input stem not specified."
	exit 1;
endif

goto check_params_return;


################################
# Error message
################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1;

arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1;


################################
# usage
################################
usage_exit:

	echo ""
	echo "USAGE: CBIG_create_subject_surf_list.csh"
	echo ""
	echo "  Required arguments"
	echo "    -sd            sub_dir      : fMRI subjects directory"
	echo "    -sub_ls        sub_list     : subjects list (only contains subjects' id, but not the full path), e.g. '/mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/surf2surf_profile_scrub_test/scripts/Subjectlist_subset'"
	echo "    -out_dir       out_dir      : output_directory"
	echo "    -surf_stem     surf_stem    : a stem that can identify the surface data that you want to use (e.g. the part that behind '*bld002', '*bld003', ...; without extension)"
	echo "  Optional arguments"
	echo "    -outlier_stem  outlier_stem : a stem that can identify the motion outliers file (without extension), e.g. '_FDRMS0.2_DVARS50_motion_outliers' (the part after '*bld002', '*bld003'; without extension)"
	echo "    -preproc_opt   preproc_opt  : assumption of preprocessing approach, choose from 'old' and 'new', 'old' means procsurffast file structure, 'new' means CBIG_fMRI_preprocess file structure. Default is 'new'"
	echo ""

	if ( $PrintHelp == 0 ) exit 1
	echo $VERSION
	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

  This function creates the surface data list and motion outliers files list. In surface data list, each line contains the full paths of one subject's surface data (only lh). The output surface list is ${sub_dir}/scripts/surf_${surf_stem}.list. In motion outliers files list, each line contains the full paths of one subject's motion outlier files. The output motion outlier files list is ${sub_dir}/scripts/outlier_${outlier_stem}.list.
  It assumes that the preprocessed surface data are located in '${sub_dir}/${subject}/surf/' and motion outliers file is located in '${sub_dir}/${subject}/qc/'. If not, the user could create symbol links to these positions.
  -surf_stem is used to distinguish the data that you want to input from other surface data. 
  -outlier_stem is similar to -surf_stem. It is used to identify the motion outliers file. For example, for surface data preprocessed by 'CBIG_fMRI_preprocess.csh', the file name of motion outliers is '*_bld*_FDRMS0.2_DVARS50_motion_outliers.txt', then outlier_stem can be '_FDRMS0.2_DVARS50_motion_outliers'. This flag is optional. If it is specified, high motion frames (outliers) will be ignored when compute functional connectivity profiles.
  
