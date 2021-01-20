#!/bin/csh -f

# example:
#     CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles.csh -lh_in <lh_profile.txt> \
#     -rh_in <rh_profile.txt> -n 7 -out <output_file>
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles.csh v 1.0 2016/06/18 $'

set lh_profile_txt = ""
set rh_profile_txt = ""
set num_clusters = ""
set output_file = ""

set mesh = fsaverage5
set mask = cortex 
set smooth = 0;
set num_tries = 1000;
set znorm = 0;
set max_iter = 100;
set no_silhouette = 1;

set PrintHelp = 0;
set n = `echo $argv | grep -e-help | wc -l`
if( $#argv == 0 || $n != 0 ) then
	echo $VERSION
	# print help	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
	exit 0;
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

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
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
	set cmd = ( $MATLAB -nojvm -nodesktop -nosplash -r '"')
	set cmd = ($cmd 'addpath(genpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"utilities"'", "'"matlab"'"')));')
	set cmd = ($cmd CBIG_AvgFreeSurferVolumes "'"${lh_profile_txt}"'" "'"${lh_avg_profile}"'"';')
	set cmd = ($cmd CBIG_AvgFreeSurferVolumes "'"${rh_profile_txt}"'" "'"${rh_avg_profile}"'"'; exit;''"')
	eval $cmd
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
	set cmd = ( $MATLAB -nojvm -nodesktop -nosplash -r '"')
	set cmd = ($cmd 'addpath(genpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"utilities"'", "'"matlab"'"')));')
	set cmd = ($cmd 'addpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"external_packages"'",)
	set cmd = ($cmd "'"SD"'", "'"SDv1.5.1-svn593"'", "'"BasicTools"'"'));')
	set cmd = ($cmd 'addpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"external_packages"'",)
	set cmd = ($cmd "'"matlab"'", "'"default_packages"'", "'"DSP"'"'));')
	set cmd = ($cmd CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun "'"${mesh}"'" "'"${mask}"'")
	set cmd = ($cmd "'"${num_clusters}"'" "'"${output_file}"'" "'"${lh_avg_profile}"'" "'"${rh_avg_profile}"'")
	set cmd = ($cmd "'"${smooth}"'" "'"${num_tries}"'" "'"${znorm}"'"' "'"${max_iter}"'"')
	set cmd = ($cmd "'"${no_silhouette}"'"'; exit;''"')
	eval $cmd
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
if( $mesh == fsaverage5 ) then
	if( $num_clusters == 17 ) then
		set ref_file = ${root_dir}/1000subjects_reference/1000subjects_clusters017_ref.mat
	else if( $num_clusters == 7 ) then
		set ref_file = ${root_dir}/1000subjects_reference/1000subjects_clusters007_ref.mat
	endif
else if( $mesh == fsaverage6 ) then
	if( $num_clusters == 17 ) then
		set ref_file = ${root_dir}/1000subjects_reference/1000subjects_clusters017_ref_NNinterp_fs6.mat
	else if( $num_clusters == 7 ) then
		set ref_file = ${root_dir}/1000subjects_reference/1000subjects_clusters007_ref_NNinterp_fs6.mat
	endif
endif

if( $mesh == fsaverage5 || $mesh == fsaverage6 ) then
	set cmd = ( $MATLAB -nojvm -nodesktop -nosplash -r '"')
	set cmd = ($cmd 'addpath(genpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"utilities"'", "'"matlab"'"')));')
	set cmd = ($cmd 'addpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"external_packages"'",)
	set cmd = ($cmd "'"matlab"'", "'"default_packages"'", "'"others"'"'));')
	set cmd = ($cmd CBIG_HungarianClusterMatchSurfWrapper "'"${ref_file}"'" "'"${output_file}"'")
	set cmd = ($cmd "'"${output_file}"'"'; exit;''"')
	eval $cmd
	echo "Hungarian match finished."
else
	echo "WARNING: The surface resolution is different from the reference file. Hungarian match is skipped."
endif


exit 0


##########################################
# parse arguments
##########################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		# surface mesh name
		case "-mesh":
			if( $#argv == 0 ) goto arg1err
			set mesh = $argv[1]; shift;
			breaksw
		
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





exit 0

#-------- Everything below is printed as part of help --------#
BEGINHELP

NAME:
	CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles.csh

DESCRIPTION:
	This function performs the group-level surface parcellation by the method of Yeo et al. 2011.
	This function will 
	(1) Average the correlation profiles of all subjects
	(2) Perform clustering algorithm on the averaged correlation profile
	(3) If the number of clusters equals to 7 or 17, the clustering results will be matched with 
	    the parcellation in Yeo et al. 2011 by Hungarian matching.
  
REQUIRED ARGUMENTS:
	-lh_in       lh_profile_txt : left hemisphere correlation profile input list (full path). 
	                              Each line in the list is the file name of the lh correlation 
	                              profile of one subject.
	-rh_in       rh_profile_txt : right hemisphere correlation profile input list (full path). 
	                              Each line in the list is the file name of the rh correlation 
	                              profile of one subject.
	-n           num_clusters   : number of clusters
	-out         output_file    : clustering output filename (full path). For example, if the 
	                              output filename is <output-dir>/cluster_007.mat, then the 
	                              user should pass in "-out <output_dir>/cluster_007". Please 
	                              be noticed that ".mat" is not included.
	
OPTIONAL ARGUMENTS:
	-tries       num_tries      : number of different random initialization (default is 1000)
	
OUTPUTS:
	<output_file>.mat
	The clustering result (.mat) file.
	
	In the same folder, there are averaged surface correlation profiles (NIFTI files):
	e.g., "lh.*.avg_profiles017.nii.gz", "rh.*.avg_profiles017.nii.gz";
	
EXAMPLE:
	CBIG_cluster_fcMRI_surf2surf_profiles.csh -lh_in \
	~/storage/fMRI_clustering/clustering_017_scrub_lh_profile.txt \
	-rh_in  ~/storage/fMRI_clustering/clustering_017_scrub_rh_profile.txt -n 17 \
	-out ~/storage/fMRI_clustering/clustering_017_scrub -tries 1000

