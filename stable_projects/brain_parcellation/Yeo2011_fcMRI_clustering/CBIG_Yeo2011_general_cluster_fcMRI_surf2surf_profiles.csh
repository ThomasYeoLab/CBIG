#!/bin/csh -f 

# example:
#     CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles.csh -sd <sub_dir> \
#     -sub_ls <sub_list.txt> -surf_stem fsaverage5 -n 7 -cluster_out GSP_100_low_motion_clusters

# This function is the main function that calls CBIG_Yeo2011_create_subject_surf_list.csh, 
# CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles_subjectlist.csh, and 
# CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles_subjectlist.csh in sequence. It assumes that the 
# preprocessed surface data are located in sub_dir/subject/surf/. surf_stem is used to distinguish 
# the data that you want to input from other surface data. 

# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

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

#########################################
# create log file
#########################################
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs
if( $scrub_flag == 0 ) then
	set LF = ${out_dir}/logs/CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles_ncluster${num_clusters}_noscrub.log
else
	set LF = ${out_dir}/logs/CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles_ncluster${num_clusters}_scrub.log
endif
rm -f $LF
touch $LF
echo "[Clustering]: logfile = $LF"


#########################################
# main function
#########################################
set cmd = "${root_dir}/CBIG_Yeo2011_create_subject_surf_list.csh -sd ${sub_dir} -sub_ls ${sub_list} -surf_stem "
set cmd = "$cmd ${surf_stem} -out_dir ${out_dir} -preproc_opt ${preproc_opt}"
if( $scrub_flag == 1 ) then
	set cmd = "$cmd -outlier_stem ${outlier_stem}"
endif
echo $cmd |& tee -a $LF
eval $cmd |& tee -a $LF

set cmd = "${root_dir}/CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list}"
set cmd = "$cmd -surf_ls ${out_dir}/lists/lh.surf${surf_stem}.list -target $target -roi $roi"
if( $scrub_flag == 1 ) then
	set cmd = "$cmd -outlier_ls ${out_dir}/lists/outlier${outlier_stem}.list"
endif
echo $cmd |& tee -a $LF
eval $cmd |& tee -a $LF

set cmd = "${root_dir}/CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list}"
set cmd = "$cmd -n ${num_clusters} -out ${cluster_out} -tries ${num_tries} -mesh $target"
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




exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

NAME:
	CBIG_Yeo2011_general_cluster_fcMRI_surf2surf_profiles.csh
	
DESCRIPTION:
	Wrapper function for group-level parcellation (Yeo et al. 2011).
	
	This function will:
	(1) Call CBIG_Yeo2011_create_subject_surf_list.csh
	    to create two surface file lists (lh and rh), where each line contains the surface files 
	    of all runs of one subject.
	(2) Call CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles_subjectlist.csh
	    To compute the correlation matrix on surface for each subject in a given subject list. 
	    When computing the correlation matrix, the users can choose to exclude the high-motion 
	    frames (motion outliers).
	(3) Call CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles_subjectlist.csh
	    To first compute the average correlation matrix across subjects. Then perform the 
	    clustering algorithm.
	
	NOTE: It assumes that the preprocessed surface data are located in '${sub_dir}/${subject}/surf/' 
	      and motion outliers file is located in '${sub_dir}/${subject}/qc/'. If not, the user could 
	      create symbol links to these positions.
	      This function will create a folder called "surf2surf_profiles" in the preprocessed fMRI 
	      folder of each subject. Therefore, write permission to the preprocessed fMRI folder of 
	      each subject is needed.
	      

REQUIRED ARGUMENTS:
	-sd            sub_dir      : fMRI subject directory. This directory contains all the folders
	                              named by the subject IDs.
	-sub_ls        sub_list     : subject list (full path). Each line in this file is one subject ID.
	-surf_stem     surf_stem    : a stem that can identify the surface data that you want to use.
	                              For example, if the surface filename is 
	                              "Sub0001_Ses1_bld002_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_
	                              sm6_fs5.nii.gz",
	                              <surf_stem> = "_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_
	                              fs5".
	-n             num_clusters : number of clusters
	-out_dir       out_dir      : output directory. It contains the log file, the lh & rh surface file 
	                              lists, and the motion outlier list. Usually, <out_dir> is the directory 
	                              pointing to <cluster_out>.
	-cluster_out   cluster_out  : clustering output filename (full path). For example, if the output 
	                              filename is <output-dir>/cluster_007.mat, then the user should pass 
	                              in "-out <output_dir>/cluster_007". Please be noticed that ".mat" 
	                              is not included. The averaged surface to surface correlation files 
	                              will be stored in the same directory as <cluster_out>.

OPTIONAL ARGUMENTS:
	-tries         num_tries    : number of different random initialization (default is 1000)
	-outlier_stem  outlier_stem : a stem that can identify the motion outliers file (without extension), 
	                              e.g. if the motion outlier file is 
	                              'Sub0001_Ses1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt',
	                              <outlier_stem> =  '_FDRMS0.2_DVARS50_motion_outliers'.
	                              If it is specified, high motion frames (outliers) will be ignored when 
	                              computing functional connectivity profiles.
	-target        target       : the resolution of clustering (default is fsaverage5)
	-roi           roi          : the resolution of ROIs (default is fsaverage3). By default, the program 
	                              will compute correlation profiles that have resolution of fsaverage5 x fsaverage3.
	-preproc_opt   preproc_opt  : assumption of preprocessing approach, choose from 'old' and 'new', 'old' 
	                              means procsurffast file structure, 'new' means CBIG_preproc_fMRI_preprocess file 
	                              structure. Default is 'new'.

OUTPUTS:
	1. A "surf2surf_profiles" folder in the preprocessed fMRI folder of each subject.
	   In this folder, 
	   
	   "*h.<subject_id>.roifsaverage3.thres0.1.surf2surf_profile.input" 
	   is a text file where each line is the surface fMRI data (lh or rh) of one run of <subject_id>.
	   
	   "outlier.<subject_id>.roifsaverage3.thres0.1.surf2surf_profile.input"
	   is a text file where each line is the motion outlier filename of one run of <subject_id>
	   
	   "*h.<subject_id>.roifsaverage3.thres0.1.surf2surf_profile_scrub.nii.gz"
	   is the surface to surface correlation file (lh or rh) of <subject_id>.
	   
	   Please see function "CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles_subjectlist.csh" for more information.
	   
	2. <out_dir>/lists/*h.surf<surf_stem>.list
	   Each line of this file is all the surface fMRI filenames (lh or rh) of one subject.
	   
	   <out_dir>/lists/outlier<outlier_stem>.list
	   Each line of this file is all the motion outlier filenames of one subject.
	   
	3. <cluster_out>.mat
	   The clustering result (.mat) file.
	   
	   In the same folder, there are averaged surface correlation profiles (NIFTI files):
	   e.g., "lh.*.avg_profiles017.nii.gz";
	   
	   and the input lists to create the averaged profiles:
	   e.g., "*_lh_profile.txt"
	   where each line is the correlation filename of one subject.

EXAMPLE:
	csh CBIG_general_cluster_fcMRI_surf2surf_profiles.csh -sd ~/storage/fMRI_data -sub_ls 
	~/storage/fMRI_data/scripts/sub_list.txt -surf_stem 
	_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5 -outlier_stem 
	_FDRMS0.2_DVARS50_motion_outliers -n 17 -out_dir ~/storage/fMRI_clustering
	-cluster_out ~/storage/fMRI_clustering/clustering_017_scrub -tries 1000 -target fsaverage5 
	-roi fsaverage3 -preproc_opt new
  
