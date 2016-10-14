#!/bin/csh -f

# Author: Jingwei Li, Date: 2016/06/18

set VERSION = '$Id: CBIG_compute_fcMRI_surf2surf_profiles.csh v 1.0 2016/06/18 $'

set sub_dir = ""
set sub = ""
set surf_data = ""
set roi = fsaverage3
set target = fsaverage5
set threshold = 0.1
set scrub_flag = 0;

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

set code_dir = `pwd`;

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print os.path.realpath('$0')"`
set root_dir = `dirname $root_dir`

set MATLAB = `which $CBIG_MATLAB_DIR/bin/matlab`
if($status) then
	echo "ERROR: could not find matlab"
	exit 1
endif

if( $scrub_flag == 0 ) then
	set output_file1 = "${output_dir}/lh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.nii.gz"
	set output_file2 = "${output_dir}/rh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.nii.gz"
else
	set output_file1 = "${output_dir}/lh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile_scrub.nii.gz"
	set output_file2 = "${output_dir}/rh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile_scrub.nii.gz"
endif

echo "===>> Compute correlation profile for ${sub_dir}/${sub}, writing outputs in $output_dir"

############################
# create input files
############################
set lh_input_file = "${output_dir}/lh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.input"
if( -e ${lh_input_file} ) then
	rm  $lh_input_file
endif
foreach surf ($surf_data)
	set lh_surf = `echo $surf | sed "s/rh/lh/g"`
	echo $lh_surf >> $lh_input_file
end

set rh_input_file = "${output_dir}/rh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.input"
if( -e ${rh_input_file} ) then
	rm  $rh_input_file
endif
foreach surf ($surf_data)
	set rh_surf = `echo $surf | sed "s/lh/rh/g"`
	echo $rh_surf >> $rh_input_file
end

if( $scrub_flag == 1 ) then
	set outlier_input_file = "${output_dir}/outlier.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.input"
	if( -e ${outlier_input_file} ) then
		rm  $outlier_input_file
	endif
	foreach outlier ($outlier_files)
		echo $outlier >> $outlier_input_file
	end
	echo "outlier_input_file = $outlier_input_file"
endif


###########################
# compute correlation profile
###########################

if( -e ${output_file1} && -e ${output_file2} ) then
	echo "Outputs already exist. Skipping......"
else
	if( $scrub_flag == 0 ) then
		$MATLAB -nodesktop -nodisplay -nosplash -r "addpath(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'surf')); CBIG_ComputeCorrelationProfile '${roi}' '${target}' '${output_file1}' '${output_file2}' '${threshold}' '${lh_input_file}' '${rh_input_file}'; exit;"
	endif
	
	if( $scrub_flag == 1 ) then
		$MATLAB -nodesktop -nodisplay -nosplash -r "addpath(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'surf')); CBIG_ComputeCorrelationProfile '${roi}' '${target}' '${output_file1}' '${output_file2}' '${threshold}' '${lh_input_file}' '${rh_input_file}' '${outlier_input_file}'; exit;"
	endif
endif

if( -e ${output_file1} && -e ${output_file2} ) then
	echo "===>> Matlab succeeds to compute correlation profiles."
else
	echo "===>> Matlab fails to compute correlation profiles."
endif

echo ""
exit 0;



############################
# parse arguments
############################
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
		
		# subject
		case "-s":
			if( $#argv == 0 ) goto arg1err;
			set sub = $argv[1]; shift;
			breaksw

		# surface data
		case "-surf_data":
			if( $#argv == 0 ) goto argerr;
			set surf_data = "$argv[1]"; shift;
			breaksw
			
		# outlier files
		case "-outlier_files":
			if( $#argv == 0 ) goto arg1err;
			set outlier_files = "$argv[1]"; shift;
			set scrub_flag = 1
			breaksw
			
		# target resolution
		case "-target":
			if( $#argv == 0 ) goto arg1err;
			set target = $argv[1]; shift;
			breaksw
			
		# ROI resolution
		case "-roi":
			if( $#argv == 0 ) goto arg1err;
			set roi = $argv[1]; shift;
			breaksw

		# output directory
		case "-output_dir":
			if( $#argv == 0 ) goto arg1err;
			set output_dir = $argv[1]; shift;
			breaksw

		default:
			echo "ERROR: Flag $flag unrecognized."
			echo $cmdline
			exit 1;
			breaksw

	endsw 
end

goto parse_args_return;


############################
# check parameters
############################
check_params:

if( $#sub_dir == 0 ) then
	echo "ERROR: subject directory not specified."
	exit 1;
endif

if( $#sub == 0 ) then
	echo "ERROR: subject id not specified."
	exit 1;
endif

if( $#surf_data == 0 ) then
	echo "ERROR: surface data not specified."
	exit 1;
endif

if( $#output_dir == 0 ) then
	echo "ERROR: output directory not specified."
	exit 1;
endif

goto check_params_return;


############################
# Error message
############################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1;

argerr:
  echo "ERROR: flag $flag requires at lease one arguments"
  exit 1;


############################
# Usage exit
############################
usage_exit:

	echo ""
	echo "USAGE: CBIG_compute_fcMRI_surf2surf_profiles.csh"
	echo ""
	echo "  Required arguments"
	echo "    -sd          sub_dir      : fMRI subejcts directory"
	echo "    -s           subject      : subject id"
	echo "    -surf_data   surf_data    : all surface data filenames of one specific subject (only left hemi, or only right hemi)"
	echo "    -output_dir  output_dir   : directory to output FC profiles file"
	echo ""
	echo "  Optional arguments"
	echo "    -outlier_ls  outlier_list : motion outliers files list"
	echo "    -target        target     : the resolution of clustering (default is fsaverage5)"
	echo "    -roi           roi        : the resolution of ROIs (defaule is fsaverage3)"
	echo ""
	
	if ( $PrintHelp == 0 ) exit 1
	echo $VERSION
	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

  This function compute the functional connectivity profiles on surface for the given subject. 
  The motion outliers file is optional. If motion outliers file is passed in, the high motion frames will be ignored when computing functional connectivity. High motion frames are indicated by '0' in motion outliers files.

