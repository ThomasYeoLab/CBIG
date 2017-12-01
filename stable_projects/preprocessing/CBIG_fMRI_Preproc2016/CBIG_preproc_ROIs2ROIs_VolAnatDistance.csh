#!/bin/csh -f

###################################################
# Compute ROIs to ROIs volumetric Euclidean distances in anatomical spaces, which is used 
# for group-level analysis, e.g. QC-RSFC correlation (the correlation between quality 
# control measure and functional connectivity) vs ROIs2ROIs Euclidean distance plot.
###################################################
# Author ##########################################
# Jingwei Li
# Oct. 13, 2017
###################################################

set VERSION = '$Id: CBIG_preproc_ROIs2ROIs_VolAnatDistance.csh, v 1.0 2017/10/13 $'


set n = `echo $argv | grep -e -help | wc -l`

# if there is no arguments or there is -help option 
if( $#argv == 0 || $n != 0 ) then
	echo $VERSION
	# print help	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
	exit 0;
endif

set n = `echo $argv | grep -e -version | wc -l`
if($n != 0) then
	echo $VERSION
	exit 0;
endif 

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

###############################
# set parameters
###############################
set default_flag = 0
if ( $default_flag == 0 ) then
	set cortical_ROIs_file = "$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/"
	set cortical_ROIs_file = "${cortical_ROIs_file}Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz"
	
	set subcortical_ROIs_file = "$CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg_182x218x182.nii.gz"
	
	set output_name = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/data/ROIs2ROIs_VolAnatDistance/"
	set output_name = "${output_name}ROIs2ROIs_dist_Schaefer2018_400Parcels_17Networks_order_19aseg_FSLMNI152_1mm.mat"
else
	set cortical_ROIs_file = $cortical_ROIs_file_in
	set subcortical_ROIs_file = $subcortical_ROIs_file_in
	set output_name = $output_name_in
endif

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

###############################
# check if matlab exists
###############################
set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
if ($status) then
	echo "ERROR: could not find matlab"
	exit 1;
endif

###############################
# Do real work! Compute distance
###############################
echo "[ROIs2ROIs VolAnatDistance]: cortical_ROIs_file = $cortical_ROIs_file"
echo "[ROIs2ROIs VolAnatDistance]: subcortical_ROIs_file = $subcortical_ROIs_file"
echo "[ROIs2ROIs VolAnatDistance]: output_name = $output_name"

set cmd = ( $MATLAB -nodesktop -nodisplay -nosplash -r '"' 'addpath(genpath('"'"${root_dir}'/utilities'"'"'))'; );
set cmd = ( $cmd CBIG_preproc_compute_ROIs2ROIs_VolAnatDistance $cortical_ROIs_file $subcortical_ROIs_file $output_name; )
set cmd = ( $cmd exit; '"' )
echo $cmd
eval $cmd

if ( -e $output_name ) then
	echo "[ROIs2ROIs VolAnatDistance]: The output file is successfully generated."
endif



exit 0


##################################
##======pass the arguments
##################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		#subject name
		case "-cortical_ROIs_file":
			if ( $#argv == 0 ) goto arg1err;
			set cortical_ROIs_file_in = $argv[1]; shift;
			breaksw	
		#path to subject's folder
		case "-subcortical_ROIs_file":
			if ( $#argv == 0 ) goto arg1err;
			set subcortical_ROIs_file_in = $argv[1]; shift;
			breaksw

		case "-output_name":
			if ( $#argv == 0 ) goto arg1err;
			set output_name_in = $argv[1]; shift;
			breaksw
			
		case "-default":
			set default_flag = 1;
			breaksw
		
		default:
			echo ERROR: Flag $flag unrecognized.
			echo $cmdline
			exit 1
			breaksw
	endsw
end
goto parse_args_return;


############################################
##======check passed parameters
############################################
check_params:
			
goto check_params_return;


#######################################################			
##======Error message		
#######################################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1

argerr:
  echo "ERROR: flag $flag requires at least one argument"
  exit 1
  
  
#######################################################
# usage sxit
#######################################################
BEGINHELP

NAME:
	CBIG_preproc_ROIs2ROIs_VolAnatDistance.csh

DESCRIPTION:
	Compute ROIs to ROIs volumetric Euclidean distance matrix, which can be used for group-level analysis, for example, 
	to plot the relationship between QC-RSFC correlation (the correlation between quality control measure and functional 
	connectivity) and ROIs2ROIs Euclidean distances.
	
	This function cannot resolve the overlap between cortical ROIs and subcortical ROIs. Pleasse make sure your input ROIs 
	are not overlapped.
	
	Make sure your "cortical_ROIs_file" and "subcortical_ROIs_file" are in the same space and with the same vox2ras matrix. 
	If both of them are passed in, the code will use the vox2ras matrix of "subcortical_ROIs_file". If the vox2ras matrices of 
	the two files are not equal, the scripts will throw an error and exit.
	
ARGUMENTS:
	-cortical_ROIs_file     cortical_ROIs_file :
	 The full name of the cortical parcellation. The default one is the 400-ROIs parcellation of Schaefer et al. 2018 in MNI 
	 1mm space:
	 $CBIG_CODE_DIR/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz.
	 If you do not want to include any cortical part, you can use 'NONE'.
	 
	-subcortical_ROIs_file  subcortical_ROIs_file :
	 The full name of the subcortical parcellation. The default one is the aseg file in MNI 1mm space:
	 $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg_182x218x182.nii.gz.
	 If you do not want to include any subcortical part, you can use 'NONE'.
	 
	-output_name            output_name :
	 The full name of the output distance matrix (matlab .mat file). The output file is a structure called "distance" which 
	 contains two fields:
	 (1) "distance.distance" is an M x M distance matrix, where M is the total number of ROIs (cortical + subcortical).
	 (2) "distance.centroids" is an M x 3 matrix contains the centroids of all ROIs.
	 The default one is 
	 "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/data/ROIs2ROIs_VolAnatDistance/
	 ROIs2ROIs_dist_Schaefer2018_400Parcels_17Networks_order_19aseg_FSLMNI152_1mm.mat"
	 
	-default :
	 If you want to use all the default values of "-cortical_ROIs_file" "-subcortical_ROIs_file", and "-output_name", 
	 you need to use this option.
	                                                
	The file in the repo: "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/data/ROIs2ROIs_VolAnatDistance/
	ROIs2ROIs_dist_Schaefer2018_400Parcels_17Networks_order_19aseg_FSLMNI152_1mm.mat" is computed by using 
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_ROIs2ROIs_VolAnatDistance.csh -default
	
EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_ROIs2ROIs_VolAnatDistance.csh 
	-cortical_ROIs_file
	"$CBIG_CODE_DIR/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz"
	-subcortical_ROIs_file "$CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg_182x218x182.nii.gz"
	-output_name 
	"$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/data/ROIs2ROIs_VolAnatDistance/
	ROIs2ROIs_dist_Schaefer2018_400Parcels_17Networks_order_19aseg_FSLMNI152_1mm.mat"

Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
	
