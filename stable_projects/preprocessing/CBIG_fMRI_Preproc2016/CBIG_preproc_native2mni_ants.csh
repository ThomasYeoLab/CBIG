#!/bin/tcsh 

# Example: 
# $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2mni_ants.csh -s Sub0001_Ses1 
#-d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS -anat_d ~/storage/sMRI_preprocess -bld '002 003' -BOLD_stem 
#_rest_stc_mc_cen_resid_lp0.08 -REG_stem _rest_stc_mc_reg -down FSL_MNI_2mm -sm 6 -sm_mask 
#${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152/SubcorticalLooseMask_MNI1mm_sm6_MNI2mm_bin0.2.nii.gz 
#-final_mask ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
#
# This function does the following steps:
#     1. Project BOLD to MNI2mm space
#     2. Smooth downsampled data (specified by -sm)
#     3. Apply final mask, if any
#
# Written by Jianxiao Wu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# BOLD: basename of each run input
# boldfolder: directory of /bold
# volfolder: store all volume results (/vol)
# bold: all runs under /bold folder

set VERSION = '$Id: CBIG_preproc_native2mni_ants.csh v 1.0 2017/08/24'

set subject = ""      # subject ID
set sub_dir = ""      # directory to subjects
set anat_dir = ""     # directory to recon-all folder
set anat_s = ""       # recon-all folder
set bold = ""         # bold number, e.g. '002 003'
set BOLD_stem = ""    # BOLD stem, e.g. _rest_stc_mc_cen_FDRMS0.2_DVARS50_resid_lp0.08
set reg_stem = ""     # registration stem, e.g. _rest_stc_mc

set nocleanup = 0;
set warp_4d = 0;
set sm = 6;
set temp_2mm = "${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz"
set temp_1mm = "${FSL_DIR}/data/standard/MNI152_T1_1mm_brain.nii.gz"

# Print help or version
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

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set currdir = `pwd`
cd $sub_dir/$subject

# create log file
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_native2mni_ants.log
if( -e ${LF} ) then
	rm $LF
endif
touch $LF
echo "[native2mni_ants]: logfile = $LF"
echo "Volumetric Projection, Downsample & Smooth" >> $LF
echo "[native2mni_ants]: CBIG_preproc_native2mni_ants.csh $cmdline"   >>$LF

set boldfolder = "$sub_dir/$subject/bold"
set volfolder = "$sub_dir/$subject/vol"
echo "[native2mni_ants]: boldfolder = $boldfolder" |& tee -a $LF

pushd $boldfolder
mkdir -p $volfolder

###################################
### register subject T1 to MNI152 1mm space
###################################
echo "======== Register $subject to MNI152 1mm space ========" |& tee -a $LF	
set mri_vol = $anat_dir/$anat_s/mri/norm.mgz
set mri_nii = $volfolder/norm.nii.gz
set warp_prefix = norm_MNI1mm
set warp = $volfolder/${warp_prefix}1Warp.nii.gz
set affine = $volfolder/${warp_prefix}0GenericAffine.mat
if(-e $warp && -e $affine) then
	echo "    [native2mni_ants]: $warp and $affine already exists." |& tee -a $LF
else
	set cmd = (mri_convert $mri_vol $mri_nii) #Convert norm.mgz to nifti format
	echo $cmd |& tee -a $LF
	eval $cmd |& tee -a $LF

	set cmd = (CBIG_antsReg_vol2vol.sh -r $temp_1mm -i $mri_nii -d $volfolder -p $warp_prefix)
	echo $cmd |& tee -a $LF
	eval $cmd |& tee -a $LF

	if(-e $warp && -e $affine) then
		echo "    [native2mni_ants]: registration to MNI152 1mm finished." |& tee -a $LF
	else
		echo "    ERROR: registration to MNI152 1mm failed." |& tee -a $LF
		exit 1;
	endif
endif
echo "======== Registration of $subject to MNI152 1mm space finished. ========" |& tee -a $LF
echo "" |& tee -a $LF

#########################################################################
### project fMRI to MNI152 2mm space & smooth
#########################################################################
foreach runfolder ($bold)
	pushd $runfolder
	echo "Run: $runfolder" |& tee -a $LF
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	echo $BOLD
	
	# if no smooth operation (sm=0), the final output is ${subject}_bld${runfolder}${BOLD_stem}_MNI2mm.nii.gz
	# if smooth (sm!=0), the final output is ${subject}_bld${runfolder}${BOLD_stem}_MNI2mm_sm${sm}.nii.gz
	set output_MNI2mm = $volfolder/${BOLD}_MNI2mm.nii.gz
	if($sm > 0) then
		set final_output = $volfolder/${BOLD}_MNI2mm_sm${sm}.nii.gz
	else
		set final_output = $output_MNI2mm
	endif

	# if final mask is applied, the final output should have _finalmask postfix
        if( $?final_mask ) then
		set final_output_mask = `basename $final_output .nii.gz`
		set final_output_mask = $volfolder/${final_output_mask}_finalmask.nii.gz
	endif
	
	# check if final output already exists
	if(-e $final_output_mask) then
		echo "[native2mni_ants]: final output: $final_output_mask already exists." |& tee -a $LF
		popd
		continue
	endif
	
	# if warp_4d is not set (default=0), then split the frames
	if($warp_4d == 0) then
		# frame_dir stores the original, projected, downsampled, and smoothed single frames
		set frame_dir = $volfolder/${BOLD}_frames
		if (-e $frame_dir) then
			echo "Warning: $frame_dir already exists"
		endif
		mkdir -p $frame_dir
	
		# split 4D image to single frame image
		set output_prefix = $frame_dir/orig_frames
		set cmd = (fslsplit $BOLD.nii.gz $output_prefix -t)
		echo $cmd |& tee -a $LF
		eval $cmd
		echo "" |& tee -a $LF
	
		set frames = `ls $frame_dir/orig_frames*nii.gz`
		set nframes = $#frames
		if($nframes == 0) then
			echo "ERROR: writing 4D volume $BOLD.nii.gz to individual frames failed" |& tee -a $LF
			exit 1;
		endif
	# otherwise, only extract the first frame (for bbregister output conversion)
	else
		set first_frame = $volfolder/orig_frame0000.nii.gz
		set cmd = (fslroi $BOLD.nii.gz $first_frame 0 1)
		echo $cmd | & tee -a $LF
		eval $cmd
		echo "" |& tee -a $LF
	endif

	###################################
	### convert bbregister output to ITK format
	###################################
	set reg_prefix = ${subject}_bld${runfolder}${reg_stem}	
	set output = $reg_prefix.txt
	if(-e $reg_prefix.txt) then
		echo "    [native2mni_ants]: $reg_prefix.txt already exists." |& tee -a $LF
	else
		if(! -e $reg_prefix.dat) then
			echo "ERROR: registration file $reg_prefix.dat not exists." |& tee -a $LF
			exit 1;
		endif

		# Convert bbregister output to .mat format
		set cmd = (tkregister2 --mov $BOLD.nii.gz --targ $mri_vol --reg $reg_prefix.dat --fslregout $reg_prefix.mat --noedit) 
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF

		# Then convert to .txt format
		if($warp_4d == 0) then
			set cmd = ($CBIG_CODE_DIR/external_packages/c3d-0.8.2-Linux-x86_64/c3d_affine_tool -ref $mri_nii)
			set cmd = ($cmd -src $frame_dir/orig_frames0000.nii.gz $reg_prefix.mat -fsl2ras -oitk $reg_prefix.txt) 
		else
			set cmd = ($CBIG_CODE_DIR/external_packages/c3d-0.8.2-Linux-x86_64/c3d_affine_tool -ref $mri_nii)
			set cmd = ($cmd -src $first_frame $reg_prefix.mat -fsl2ras -oitk $reg_prefix.txt)
		endif
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF

		if($warp_4d == 1) then
			rm $first_frame
		endif

		if(! -e $reg_prefix.txt) then
			echo "ERROR: converting bbregister warp to $reg_prefix.txt failed." |& tee -a $LF
			exit 1;
		endif
	endif

	###################################
	### project fMRI to MNI152 2mm space (using 1mm warp)
	###################################
	echo "======== Project $runfolder to MNI152 2mm space ========" |& tee -a $LF
	if($warp_4d == 0) then
		# Project each frame separately
		set fcount = 0;
		while($fcount < $nframes)
			set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`		
			set input = $frame_dir/orig_frames${fcount_str}.nii.gz
			set output = $frame_dir/${fcount_str}_MNI2mm.nii.gz
			if(-e $output) then
				echo "    [native2mni_ants]: $output already exists." |& tee -a $LF
			else
				set cmd = (CBIG_antsApplyReg_vol2vol.sh -i $input -r $temp_2mm -w $warp_prefix -p)
				set cmd = ($cmd ${fcount_str}_MNI2mm -f $reg_prefix.txt -d $volfolder -o $frame_dir) 
				echo $cmd |& tee -a $LF
				eval $cmd |& tee -a $LF

				if(-e $output) then
					echo "    [native2mni_ants]: projection to $output finished." |& tee -a $LF
				else
					echo "    ERROR: projection to $output failed." |& tee -a $LF
					exit 1;
				endif
			endif
		
			@ fcount = $fcount + 1
		end
	else
		# If warp_4d is set to 1, project the 4D volume directly
		set output = $output_MNI2mm
		if (-e $output) then
			echo "    [native2mni_ants]: $output already exists." |& tee -a $LF
		else
			set cmd = (CBIG_antsApplyReg_vol2vol.sh -i $BOLD.nii.gz -r $temp_2mm -w $warp_prefix -p)
			set cmd = ($cmd ${BOLD}_MNI2mm -f $reg_prefix.txt -e 1 -d $volfolder -o $volfolder)
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
	
			if(-e $output) then
				echo "    [native2mni_ants]: projection to $output finished." |& tee -a $LF
			else
				echo "    ERROR: projection to $output failed." |& tee -a $LF
				exit 1;
			endif
		endif
	endif
	echo "======== Projection of $runfolder to MNI152 2mm space finished. ========" |& tee -a $LF
	echo "" |& tee -a $LF
	
	#########################
	### smooth
	#########################
	echo "======== Smooth in MNI152 2mm space with fwhm = $sm ========" |& tee -a $LF
	if($warp_4d == 0) then
		# Smooth each frame separately
		set fcount = 0;
		while($fcount < $nframes)
			set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`	
			set input = $frame_dir/${fcount_str}_MNI2mm.nii.gz
			mkdir -p $frame_dir/sm
			set output = $frame_dir/sm/${fcount_str}_MNI2mm_sm${sm}.nii.gz
			if(-e $output) then
				echo "[native2mni_ants]: $output already exists." |& tee -a $LF
			else
				set tmp1 = $frame_dir/tmp1_${fcount_str}.nii.gz
				set tmp2 = $frame_dir/tmp2_${fcount_str}.nii.gz
				set std = `awk "BEGIN {print ${sm}/2.35482}"`    #Note that fwhm = 2.35482 * std, fslmaths -s 
										 #is in unit of mm, not voxel.
				if($?sm_mask) then
					# if the user passes in a volume sm_mask, the procedure is
					# 1. smooth volume data within sm_mask
					# 2. smooth sm_mask within sm_mask
					# 3. divide smoothed volume by smoothed sm_mask (deal with boundary problem)
				
					set cmd = (fslmaths $input -s $std -mas $sm_mask $tmp1)
					echo $cmd |& tee -a $LF
					eval $cmd
				
					set cmd = (fslmaths $sm_mask -s $std -mas $sm_mask $tmp2)
					echo $cmd |& tee -a $LF
					eval $cmd
				
					set cmd = (fslmaths $tmp1 -div $tmp2 $output)
					echo $cmd |& tee -a $LF
					eval $cmd
				else
					set cmd = (fslmaths $input -s $std $output)
					echo $cmd |& tee -a $LF
					eval $cmd
				endif
			
				if(-e $output) then
					echo "[native2mni_ants]: smooth to $output finished." |& tee -a $LF
				else
					echo "ERROR: smooth to $output failed." |& tee -a $LF	
					exit 1;
				endif
			endif

			@ fcount = $fcount + 1
		end
	else
		# If warp_4d is set to 1, smooth the 4D volume directly
		set input = $output_MNI2mm
		set output = $final_output
		if (-e $output) then
			echo "[native2mni_ants]: $output already exists." |& tee -a $LF
		else
			set tmp1 = $volfolder/tmp1_${BOLD}_MNI2mm.nii.gz
			set tmp2 = $volfolder/tmp2_${BOLD}_MNI2mm.nii.gz
			set std = `awk "BEGIN {print ${sm}/2.35482}"`    #Note that fwhm = 2.35482 * std, fslmaths -s is in 
									 #unit of mm, not voxel.
			if($?sm_mask) then
				# if the user passes in a volume sm_mask, the procedure is
				# 1. smooth volume data within sm_mask
				# 2. smooth sm_mask within sm_mask
				# 3. divide smoothed volume by smoothed sm_mask (deal with boundary problem)
				
				set cmd = (fslmaths $input -s $std -mas $sm_mask $tmp1)
				echo $cmd |& tee -a $LF
				eval $cmd
				
				set cmd = (fslmaths $sm_mask -s $std -mas $sm_mask $tmp2)
				echo $cmd |& tee -a $LF
				eval $cmd
				
				set cmd = (fslmaths $tmp1 -div $tmp2 $output)
				echo $cmd |& tee -a $LF
				eval $cmd
			else
				set cmd = (fslmaths $input -s $std $output)
				echo $cmd |& tee -a $LF
				eval $cmd
			endif

			if( $nocleanup == 0 ) then
				rm ${input}
				rm ${tmp1}
				rm ${tmp2}
			endif
		
			if (-e $output) then
				echo "[native2mni_ants]: smooth to $output finished." |& tee -a $LF
			else
				echo "ERROR: smooth to $output failed." |& tee -a $LF	
				exit 1;
			endif
		endif
	endif
	echo "" |& tee -a $LF
	
	######################
	### combine frames if necessary
	######################
	if($warp_4d == 0) then
		echo "======== Combine frames for $runfolder ========" |& tee -a $LF
		if($sm > 0) then
			set cmd = (fslmerge -t $final_output $frame_dir/sm/*_MNI2mm_sm${sm}.nii.gz)
		else
			set cmd = (fslmerge -t $final_output $frame_dir/*_MNI2mm.nii.gz)
		endif
		echo $cmd |& tee -a $LF
		eval $cmd
	
		# for some reason, the split-merge procedure will change TR in header information
		# now change it back
		set MRI_info = `mri_info $BOLD.nii.gz`
		set TR = `echo $MRI_info | grep -o 'TR: \(.*\)' | awk -F " " '{print $2}'`
	
		set cmd = (mri_convert $final_output $final_output -tr $TR)
		echo $cmd |& tee -a $LF
		eval $cmd
	
		if(-e $final_output) then
			echo "======== Combination of $runfolder finished ========" |& tee -a $LF
		else
			echo "======== Combination of $runfolder failed ========" |& tee -a $LF
			exit 1;
		endif
		echo "" |& tee -a $LF
	endif

	#########################
	### apply final mask
	#########################
        if( $?final_mask ) then
		echo "======== Applying final mask for $runfolder ========" |& tee -a $LF
		set input = $final_output
		set output = $final_output_mask
		
		if( -e $output ) then  
  			echo "[native2mni_ants]: $output already exists" |& tee -a $LF
 		else		
			set cmd = (fslmaths $input -mas ${final_mask} $output)
			echo $cmd
			eval $cmd

			if( $nocleanup == 0 ) then
				rm ${input}
			endif
			
			if( -e $output ) then
				echo "[native2mni_ants]: Applying final mask finished. The output is $output" |& tee -a $LF
			else
				echo "ERROR: Applying final mask failed." |& tee -a $LF
			endif			
		endif
        endif
	
	
	#########################
	### clean up
	#########################
	if( $nocleanup == 0 & $warp_4d == 0) then
		rm -R $frame_dir
		echo "======== Cleaning up ========" |& tee -a $LF
		echo "[native2mni_ants]: cleaned up intermediate results in $frame_dir" |& tee -a $LF
	endif	
	
	popd
end

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2mni_ants.csh >> $LF
endif

echo "****************************************************************" |& tee -a $LF
echo "" |& tee -a $LF

popd
exit 0;

#######################################
##======pass the arguments
#######################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		# subject name (required)
		case "-s":
			if ( $#argv == 0 ) goto arg1err;
			set subject = $argv[1]; shift;
			breaksw	
		
		# path to subject's folder (required)
		case "-d":
			if ( $#argv == 0 ) goto arg1err;
			set sub_dir = $argv[1]; shift;
			breaksw
		
		# anatomical directory (required)
		case "-anat_d":
			if ( $#argv == 0 ) goto arg1err;
			set anat_dir = $argv[1]; shift;
			breaksw
		
		# anatomical data (required)
		case "-anat_s":
			if ( $#argv == 0 ) goto arg1err;
			set anat_s = $argv[1]; shift;
			breaksw

		# bold number, like "002 003"	
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set bold = ($argv[1]); shift;
			breaksw
			
		# MNI space smooth param
		case "-sm":
			if ( $#argv == 0 ) goto arg1err;
			set sm = $argv[1]; shift;
			breaksw
			
		# smooth mask
		case "-sm_mask":
			if ( $#argv == 0 ) goto arg1err;
			set sm_mask = $argv[1]; shift;
			breaksw
		
		# BOLDbase_suffix 
		case "-BOLD_stem":
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = "$argv[1]"; shift;
			breaksw
			
		# reg suffix (default rest_skip*_stc_mc_reg)
		case "-REG_stem":
			if ( $#argv == 0 ) goto arg1err;
			set reg_stem = "$argv[1]"; shift;
			breaksw
			
		# final mask (a loose MNI brain mask)
		case "-final_mask":
			if ( $#argv == 0 ) goto arg1err;
			set final_mask = $argv[1]; shift;
			breaksw
			
		# if nocleanup is turned on, it won't delete the folders with single frames, otherwise they will be removed
		case "-nocleanup":
			set nocleanup = 1
			breaksw
		
		# update results, if exist then do not generate
		case "-warp_4d":
			set warp_4d = 1;
			breaksw
		
		default:
			echo ERROR: Flag $flag unrecognized.
			echo $cmdline
			exit 1
			breaksw
	endsw
end
goto parse_args_return;


#############################################
##======check passed parameters
#############################################
check_params:
if ( $#subject == 0 ) then
	echo "ERROR: subject not specified"
	exit 1;
endif
 
if ( $#sub_dir == 0 ) then
	echo "ERROR: path to subject folder not specified"
	exit 1;
endif	

if ( $#anat_dir == 0 ) then
	echo "ERROR: anatomical directory not specified"
	exit 1;
endif

if ( $#anat_s == 0 ) then
	echo "ERROR: anatomical data not specified"
	exit 1;
endif

if ( $#bold == 0 ) then
	echo "ERROR: bold number not specified"
	exit 1;
endif

if ( $#BOLD_stem == 0 ) then
	echo "ERROR: BOLD stem not specified"
	exit 1;
endif

if ( $#reg_stem == 0 ) then
	echo "ERROR: reg stem not specified"
	exit 1;
endif
				
goto check_params_return;

#########################################			
##======Error message		
#########################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1

argerr:
  echo "ERROR: flag $flag requires at least one argument"
  exit 1

#########################################
# Usage exit
#########################################
BEGINHELP

NAME:
	CBIG_preproc_native2mni_ants.csh
	
DESCRIPTION:
	Project data from fMRI native space to MNI152 2mm space using ANTs and smooth.
	
	This function does the following steps:
	    1. Project fMRI volumes to MNI152 2mm space usig ANTs
	    2. Smooth the downsampled data. The user can specify a mask for smoothing.
	    3. Apply the final mask, if any.

	Registration step uses antsRegistration to register T1 volume to MNI152 1mm space. Then, projection 
	step uses antsApplyWarps to project fMRI data to MNI152 2mm space. 

	For smoothing step, the users can pass in a mask (using -sm_mask flag). If the users specify a mask, 
	then the program will smooth the volume within the mask, smooth the mask, and divide the smoothed 
	volume by the smoothed mask (to deal with boundary problems).
	
	Caution: check your ANTs software version before using. There is a bug in early builds of ANTs 
	(before Aug 2014) that causes resampling for timeseries to be wrong. We have tested that our codes 
	would work on ANTs version 2.2.0.

REQUIRED ARGUMENTS:
	-s          subject   : fMRI subject id
	-d          sub_dir   : absolute path to <subject>, all preprocessed data of this 
	                        subject is assumed to be stored in <sub_dir>/<subject>
	-anat_s     anat_s    : anatomical subject id (recon-all output subject name)
	-anat_d     anat_dir  : absolute path to <anat_s>. All recon-all outputs of this 
	                        subject is stored in <anat_dir>/<anat_s>.
	-bld        bold      : all bold run numbers of this subject. Each number must be 
	                        three digits. If this subject has multiple runs, use a space 
	                        as delimiter, e.g. '002 003'. NOTE: quote sign is necessary.
	-BOLD_stem  BOLD_stem : stem of input file, e.g. if the input file name is
	                        Sub0001_Ses1_bld002_rest_stc_mc_cen_resid_cen_FDRMS0.2_DVARS50_bp_0_0.08.nii.gz,
	                        then <BOLD_stem> = _rest_stc_mc_cen_resid_cen_FDRMS0.2_DVARS50_bp_0_0.08. 
	                        This input file is assumed to be in <sub_dir>/<subject>/bold/<run_number>.
	-REG_stem   reg_stem  : stem of T1-T2* registration file. E.g. if the registration 
	                        file is Sub0001_Ses1_bld002_rest_skip4_stc_mc_reg.dat, then 
	                        <REG_stem> = _rest_skip4_stc_mc_reg. The registration file 
	                        is assumed to be stored in <sub_dir>/<subject>/bold/<run_number>.
	                        
OPTIONAL ARGUMENTS:
	-sm         sm        : smooth fwhm in mm (default: 6)
	-sm_mask    sm_mask   : mask for smoothing (e.g. a grey matter mask in MNI152 2mm)
	-final_mask final_mask: define the loose mask in MNI152 2mm space applied to the volume in the
	                        last step. An example of the final mask is: 
				${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
	-nocleanup            : do not remove intermediate results.
	-warp_4d              : use this flag to project the 4D volume of frames without splitting. 
				otherwise, the frames are split and projected separately by default.
	
OUTPUTS:	
	(1) The fMRI volume after projection (in MNI152 2mm space):
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI2mm.nii.gz
	
	(2) The fMRI volume after smooth:
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI2mm_sm<sm>.nii.gz
	
	(3) The fMRI volume after applying the final mask (if -final_mask is used):
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI2mm_sm<sm>_finalmask.nii.gz
	
	(4) Intermediate results (if -nocleanup is passed in):
	    Folders <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem> including all intermediate 
	    fMRI single frames. (The projection, downsample, and smoothing are done frame by frame to save 
	    memory, if -warp_4d is not used)

Example:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2mni.csh -s Sub0001_Ses1 
	-d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS -anat_d ~/storage/sMRI_preprocess -bld '002 003' 
	-BOLD_stem _rest_stc_mc_cen_resid_lp0.08 -REG_stem _rest_stc_mc_reg -sm 6 -sm_mask 
	${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152/SubcorticalLooseMask_MNI1mm_sm6_MNI2mm_bin0.2.nii.gz 
	-final_mask ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz


