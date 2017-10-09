#!/bin/csh -f 

# Author: Jingwei Li, Date: 2016/06/18
# example:
#     CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh -sd /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast -sub_ls /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast/scripts/GSP_newtestlist.txt -n 7 -out /mnt/eql/yeo2/CBIG_repo_tests/CBIG_fMRI_Preproc2016_test/GSP_surface_100sub/procsurffast/clustering/GSP_100_low_motion_clusters

set VERSION = '$Id: CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh v 1.0 2016/06/18 $'

set sub_dir = ""
set subjects = ""
set output_file = ""
set roi = fsaverage3
set threshold = 0.1
set num_tries = 1000
set scrub_flag = 0

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

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

set output_dir = `dirname $output_file`
mkdir -p $output_dir

set lh_profile_input = ${output_file}_lh_profile.txt
rm $lh_profile_input

set rh_profile_input = ${output_file}_rh_profile.txt
rm $rh_profile_input

set subjects = `cat ${sub_list}`
foreach s ($subjects)
	if( $scrub_flag == 0 ) then
		echo "${sub_dir}/$s/surf2surf_profiles/lh.$s.roi${roi}.thres${threshold}.surf2surf_profile.nii.gz" >> $lh_profile_input
		echo "${sub_dir}/$s/surf2surf_profiles/rh.$s.roi${roi}.thres${threshold}.surf2surf_profile.nii.gz" >> $rh_profile_input
	else
		echo "${sub_dir}/$s/surf2surf_profiles/lh.$s.roi${roi}.thres${threshold}.surf2surf_profile_scrub.nii.gz" >> $lh_profile_input
		echo "${sub_dir}/$s/surf2surf_profiles/rh.$s.roi${roi}.thres${threshold}.surf2surf_profile_scrub.nii.gz" >> $rh_profile_input
	endif
end

set cmd = "${root_dir}/CBIG_cluster_fcMRI_surf2surf_profiles.csh -lh_in ${lh_profile_input} -rh_in ${rh_profile_input} -n ${num_clusters} -out ${output_file} -tries ${num_tries}"
echo $cmd
eval $cmd

exit 0

#############################
# parse arguments
#############################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		# subjects directory
		case "-sd":
			if( $#argv == 0 ) goto arg1err;
			set sub_dir = $argv[1]; shift;
			breaksw
			
		# subject list
		case "-sub_ls":
			if( $#argv == 0 ) goto arg1err;
			set sub_list = $argv[1]; shift;
			breaksw
			
		# number of clusters
		case "-n":
			if( $#argv == 0 ) goto arg1err;
			set num_clusters = $argv[1]; shift;
			breaksw
			
		# output file
		case "-out":
			if( $#argv == 0 ) goto arg1err;
			set output_file = $argv[1]; shift;
			breaksw
			
		case "-tries":
			if( $#argv == 0 ) goto arg1err
			set num_tries = $argv[1]; shift;
			breaksw
			
		case "-scrub_flag":
			if( $#argv == 0 ) goto arg1err
			set scrub_flag = $argv[1]; shift;
			breaksw
			
		default:
			echo "ERROR: flag $flag unrecognized."
			echo $cmdline
			exit 1
			breaksw
		
	endsw
end

goto parse_args_return;


#############################
# check parameters
#############################
check_params:

if( $#sub_dir == 0 ) then
	echo "ERROR: subjects directory not specified."
	exit 1;
endif

if( $#sub_list == 0 ) then
	echo "ERROR: subject list not specified."
	exit 1;
endif

if( $#num_clusters == 0 ) then
	echo "ERROR: number of clusters not specified."
	exit 1;
endif

if( $#output_file == 0 ) then
	echo "ERROR: output file not specified."
	exit 1;
endif

if( $scrub_flag != 0  && $scrub_flag != 1 ) then
	echo "ERROR: wrong input for scrub_flag."
	exit 1;
endif

goto check_params_return;


##############################
# Error message
##############################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
  



exit 0

#-------- Everything below is printed as part of help --------#
BEGINHELP

NAME:
	CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh

DESCRIPTION:
	This function is the wrapper function for clustering. It generates the lists that contain 
	the surface profiles filenames of all subjects. And it calls 
	"CBIG_cluster_fcMRI_surf2surf_profiles.csh" to perform clustering.
  
REQUIRED ARGUMENTS:
	-sd          sub_dir      : fMRI subjects directory. This directory contains all the folders
	                            named by the subject IDs.
	-sub_ls      sub_list     : subjects list (full path). Each line in this file is one subject ID.
	-n           num_clusters : number of clusters
	-out         output_file  : clustering output filename (full path). For example, if the output 
	                            filename is <output-dir>/cluster_007.mat, then the user should pass 
	                            in "-out <output_dir>/cluster_007". Please be noticed that ".mat" 
	                            is not included.

OPTIONAL ARGUMENTS:
	-tries       num_tries    : number of different random initialization for clustering (default is 1000)
	-scrub_flag  scrub_flag   : 0 or 1, 1 for ignoring high motion frames when computing profiles; 
	                            0 for keeping all frames when computing correlation profiles. 
	                            Default is 0. This option is used to differentiate the filenames 
	                            for outputs.

OUTPUTS:
	<output_file>.mat
	The clustering result (.mat) file.
	
	In the same folder, there are averaged surface correlation profiles (NIFTI files):
	e.g., "lh.*.avg_profiles017.nii.gz";
	
	and the input lists to create the averaged profiles:
	e.g., "*_lh_profile.txt"
	where each line is the correlation filename of one subject.
	
EXAMPLE:
	csh CBIG_cluster_fcMRI_surf2surf_profiles_subjectlist.csh -sd ~/storage/fMRI_data -sub_ls 
	~/storage/fMRI_data/scripts/sub_list.txt -n 17 -out ~/storage/fMRI_clustering/clustering_017_scrub 
	-tries 1000 -scrub_flag 1


