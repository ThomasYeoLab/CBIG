#!/bin/csh -f

# Author : Jingwei Li, Date: 2016/06/18
# example:
#     CBIG_cluster_fcMRI_surf2surf_profiles.csh -lh_in /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast/clustering/GSP_100_low_motion_clusters_lh_profile.txt -rh_in /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast/clustering/GSP_100_low_motion_clusters_rh_profile.txt -n 7 -out /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast/clustering/GSP_100_low_motion_clusters

set VERSION = '$Id: CBIG_cluster_fcMRI_surf2surf_profiles.csh v 1.0 2016/06/18 $'

set lh_profile_txt = ""
set rh_profile_txt = ""
set num_clusters = ""
set output_file = ""

set mesh = fsaverage5
set mask = cortex 
set smooth = 0;
set num_tries = 1000;
set znorm = 0;

set PrintHelp = 0;
if( $#argv == 0 ) goto usage_exit;
set n = `echo $argv | grep -e-help | wc -l`
if( $n != 0 ) then
	set PrintHelp = 1;
	goto usage_exit;
endif
set n = `echo $argv | grep -e-version | wc -l`
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

set MATLAB = `which $CBIG_MATLAB_DIR/bin/matlab`
if($status) then
	echo "ERROR: could not find matlab"
	exit 1
endif

set output_dir = `dirname $output_file`
echo "output_dir = $output_dir"

#######################################
# Average profiles
#######################################
echo "===>> Average profiles"
set formated_cluster = `echo $num_clusters | awk '{printf ("%03d", $1)}'`
set output_base = `basename ${output_file}`
set lh_avg_profile = ${output_dir}/lh.${output_base}.avg_profiles${formated_cluster}.nii.gz
set rh_avg_profile = ${output_dir}/rh.${output_base}.avg_profiles${formated_cluster}.nii.gz

echo "lh_avg_profile = $lh_avg_profile"
echo "rh_avg_profile = $rh_avg_profile"

if( -e $lh_avg_profile && -e $rh_avg_profile ) then
	echo "Averaged profiles $lh_avg_profile & $rh_avg_profile already exist. Skipping ......"
else
	$MATLAB -nojvm -nodesktop -nosplash -r "addpath(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'utilities')); CBIG_AvgFreeSurferVolumes '${lh_profile_txt}' '${lh_avg_profile}'; CBIG_AvgFreeSurferVolumes '${rh_profile_txt}' '${rh_avg_profile}'; exit;"
endif

if( ! -e $lh_avg_profile ) then
	echo "ERROR: $lh_avg_profile not produced."
	exit 1
endif
if( ! -e $rh_avg_profile ) then
	echo "ERROR: $rh_avg_profile not produced."
	exit 1
endif

echo "===>> Averaging profiles finished."
echo ""


########################################
# clustering
########################################
echo "===>> Cluster"

if( -e "${output_file}.mat" ) then
	echo "Clustering results ${output_file}.mat already exist. Skipping ......"
else
	$MATLAB -nojvm -nodesktop -nosplash -r "addpath(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'DSP')); CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun '$mesh' '$mask' '${num_clusters}' '${output_file}' '${lh_avg_profile}' '${rh_avg_profile}' '$smooth' '${num_tries}' '${znorm}'; exit;"
endif

if( ! -e "${output_file}.mat" ) then
	echo "ERROR: clustering unsuccesful. ${output_file}.mat not produced."
	exit 1
endif

echo "===>> Clustering finished."
echo ""


########################################
# Hungarian match
########################################
echo "===>> Hungarian match"
if( $num_clusters == 17 ) then
	set ref_file = ${root_dir}/1000subjects_reference/1000subjects_clusters017_ref.mat
else if( $num_clusters == 7 ) then
	set ref_file = ${root_dir}/1000subjects_reference/1000subjects_clusters007_ref.mat
endif

$MATLAB -nojvm -nodesktop -nosplash -r "addpath(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'utilities')); CBIG_HungarianClusterMatchSurfWrapper '${ref_file}' '${output_file}' '${output_file}'; exit;"

echo "Hungarian match finished."

exit 0


##########################################
# parse arguments
##########################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		# lh input
		case "-lh_in":
			if( $#argv == 0 ) goto arg1err
			set lh_profile_txt = $argv[1]; shift;
			breaksw
			
		# rh input
		case "-rh_in":
			if( $#argv == 0 ) goto arg1err
			set rh_profile_txt = $argv[1]; shift;
			breaksw
			
		# number of clusters
		case "-n":
			if( $#argv == 0 ) goto arg1err
			set num_clusters = $argv[1]; shift;
			breaksw
		
		# output file
		case "-out":
			if( $#argv == 0 ) goto arg1err
			set output_file = $argv[1]; shift;
			breaksw
			
		case "-tries":
			if( $#argv == 0 ) goto arg1err
			set num_tries = $argv[1]; shift;
			breaksw
			
		default:
			echo "ERROR: flag $flag unrecognized."
			echo $cmdline
			exit 1
			breaksw
		
	endsw
end

goto parse_args_return;


##########################################
# check parameters
##########################################
check_params:

if( $#lh_profile_txt == 0 ) then
	echo "ERROR: left hemisphere profile input not specified."
	exit 1
endif

if( $#rh_profile_txt == 0 ) then
	echo "ERROR: right hemisphere profile input not specified."
	exit 1
endif

if( $#num_clusters == 0 ) then
	echo "ERROR: number of clusters not specified."
	exit 1
endif

if( $#output_file == 0 ) then
	echo "ERROR: output file not specified."
	exit 1
endif

goto check_params_return;


##########################################
# Error message
##########################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1


##########################################
# Usage exit
##########################################
usage_exit:

	echo ""
	echo "USAGE: CBIG_cluster_fcMRI_surf2surf_profiles.csh"
	echo ""
	echo "  Required arguments"
	echo "    -lh_in       lh_profile_txt : left hemisphere profile input"
	echo "    -rh_in       rh_profile_txt : right hemisphere profile input"
	echo "    -n           num_clusters   : number of clusters"
	echo "    -out         output_file    : clustering output filename (full path, without extension)"
	echo ""
	echo "  Optional arguments"
	echo "    -tries       num_tries      : number of different random initialization (default is 1000)"
	echo ""

	if ( $PrintHelp == 0 ) exit 1
	echo $VERSION
	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

  This function is called by 'CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh'. It needs to take in left and right hemispheres FC profiles lists. The two lists are generated inside 'CBIG_cluster_fcMRI_surf2surf_profiles.csh'. This is the main function for clustering.
  
