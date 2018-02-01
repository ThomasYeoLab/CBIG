#!/bin/csh -f

# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$ID: CBIG_Yeo2011_create_subject_surf_list.csh, v 1.0 2016/06/18'

set sub_dir = "";
set sub_list = "";
set fsaverage = "";
set surf_stem = ""
set outlier_stem = ""
set scrub_flag = 0;
set out_dir = ""
set preproc_opt = "new"

set PrintHelp = 0;
set n = `echo $argv | grep -e -help | wc -l`
if( $#argv == 0 || $n != 0 ) then
	echo $VERSION
	# print help	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
	exit 0;
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


	

exit 0

#-------- Everything below is printed as part of help --------#
BEGINHELP

NAME:
	CBIG_Yeo2011_create_subject_surf_list.csh

DESCRIPTION:
	This function creates the surface fMRI data list and motion outliers files list. 
	In surface data list, each line contains the full paths of one subject's surface data (lh or rh). 
	In motion outliers files list, each line contains the full paths of one subject's motion outlier files. 
	
	It assumes that the preprocessed surface data are located in '${sub_dir}/${subject_id}/surf/' and the
	motion outliers file is located in '${sub_dir}/${subject_id}/qc/'. If not, the user could create symbol 
	links to these locations.

REQUIRED ARGUMENTS:
	-sd            sub_dir      : fMRI subjects directory. This directory contains all the folders
	                              named by the subject IDs.
	-sub_ls        sub_list     : subject list (full path). Each line in this file is one subject ID.
	-out_dir       out_dir      : output directory. It contains the log file, the lh & rh surface fMRI 
	                              file lists, and the motion outlier list. 
	-surf_stem     surf_stem    : a stem that can identify the surface data that you want to use.
	                              For example, if the surface file name is 
	                              "Sub0001_Ses1_bld002_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz",
	                              <surf_stem> = "_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5".
	
OPTIONAL ARGUMENTS:
	-outlier_stem  outlier_stem : a stem that can identify the motion outliers file (without extension), 
	                              e.g. if the motion outlier file is 
	                              'Sub0001_Ses1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt',
	                              <outlier_stem> = '_FDRMS0.2_DVARS50_motion_outliers'.
	                              If it is not specified, the motion outlier file name list will not be created.
	-preproc_opt   preproc_opt  : assumption of preprocessing approach, choose from 'old' and 'new', 'old' 
	                              means procsurffast file structure, 'new' means CBIG_preproc_fMRI_preprocess file 
	                              structure. Default is 'new'.

OUTPUTS:
	<out_dir>/lists/lh.surf<surf_stem>.list
	Each line of this file is all the surface fMRI file names (lh) of one subject.
	
	<out_dir>/lists/rh.surf<surf_stem>.list
	Each line of this file is all the surface fMRI file names (rh) of one subject.
	
	<out_dir>/lists/outlier<outlier_stem>.list
	Each line of this file is all the motion outlier file names of one subject.
	
EXAMPLE:
	csh CBIG_create_subject_surf_list.csh -sd ~/storage/fMRI_data -sub_ls 
	~/storage/fMRI_data/scripts/sub_list.txt -surf_stem 
	_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5 -outlier_stem 
	_FDRMS0.2_DVARS50_motion_outliers -out_dir ~/storage/fMRI_clustering  -preproc_opt new
  
