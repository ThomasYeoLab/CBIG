#!/bin/tcsh 
 
# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2mni.csh 
#	-s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS -anat_d ~/storage/sMRI_preprocess 
#	-bld '002 003' -BOLD_stem _rest_stc_mc_cen_resid_lp0.08 -REG_stem _rest_stc_mc_reg -down FSL_MNI_2mm -sm 6 -sm_mask
#	${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152_masks/SubcorticalLooseMask_MNI1mm_sm6_MNI2mm_bin0.2.nii.gz
#	-final_mask ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
#
# This function does the following steps:
#     1. Project anatomical data to MNI1mm space for check purpose
#     2. Project BOLD to MNI1mm space
#     3. Downsample projected data to MNI2mm (either FS 2mm or FSL 2mm, specified by -down)
#     4. Smooth downsampled data (specified by -sm)
#     5. Apply final mask, if any
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#BOLD: basename of each run input
#boldfolder: directory of /bold
#volfolder: store all volume results (/vol)
#bold: all runs under /bold folder

set VERSION = '$Id: CBIG_preproc_native2mni.csh v 1.0 2016/05/30'

set subject = ""      # subject ID
set sub_dir = ""      # directory to subjects
set anat_dir = ""     # directory to recon-all folder
set anat_s = ""       # recon-all folder
set bold = ""         # bold number, e.g. '002 003'
set BOLD_stem = ""    # BOLD stem, e.g. _rest_stc_mc_cen_FDRMS0.2_DVARS50_resid_lp0.08
set reg_stem = ""     # registration stem, e.g. _rest_stc_mc
set down = ""         # downsample flag, FSL_MNI_FS_2mm for FreeSurfer 2mm space, FSL_MNI_2mm for FSL MNI 2mm space

set force = 0;
set nocleanup = 0;
set sm = 7;
set FS_temp_2mm = "${CBIG_CODE_DIR}/data/templates/volume/FS_nonlinear_volumetric_space_4.5/gca_mean2mm.nii.gz"
set MNI_temp_2mm = "${FSL_DIR}/data/standard/MNI152_T1_2mm_brain.nii.gz"
set temp_2mm = ${MNI_temp_2mm}
set MNI_ref_dir = "${CBIG_CODE_DIR}/data/templates/volume"
set MNI_ref_id = "FSL_MNI152_FS4.5.0"

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


setenv SUBJECTS_DIR $anat_dir

set currdir = `pwd`
cd $sub_dir/$subject


# create log file
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_native2mni.log
if( -e ${LF} ) then
	rm $LF
endif
touch $LF
echo "[native2mni]: logfile = $LF"
echo "Volumetric Projection, Downsample & Smooth" >> $LF
echo "[native2mni]: CBIG_preproc_native2mni.csh $cmdline"   >>$LF

set boldfolder = "$sub_dir/$subject/bold"
set volfolder = "$sub_dir/$subject/vol"
echo "[native2mni]: boldfolder = $boldfolder" |& tee -a $LF

pushd $boldfolder
mkdir -p $volfolder

set fs_version = `cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' | head -c 1`
echo "Freesurfer version: ${fs_version}" |& tee -a $LF


#############################################################
### project T1 to MNI152 1mm space for check purpose
#############################################################
echo "================== Project T1 to MNI152 1mm space ==================" |& tee -a $LF
set input = ${anat_dir}/${anat_s}/mri/norm.mgz
set output = $volfolder/norm_MNI152_1mm.nii.gz
if(-e $output) then
	echo "[native2mni]: $output already exists." |& tee -a $LF
else
	set cmd = (CBIG_vol2vol_m3z.csh -src-id $anat_s -src-dir $anat_dir -targ-id $MNI_ref_id -targ-dir $MNI_ref_dir -in $input -out $output -no-cleanup)
	echo $cmd |& tee -a $LF
	$cmd |& tee -a $LF
	if(-e $output) then
		echo "========== Projection of T1 to MNI152 1mm space finished ==========" |& tee -a $LF
	else
		echo "ERROR: Projection T1 to MNI152 1mm space failed." |& tee -a $LF
		exit 1;
	endif
endif
if( $nocleanup == 0 ) then
	rm $volfolder/FStmp.norm_MNI152_1mm.nii.gz
endif
echo "" |& tee -a $LF


#########################################################################
### project fMRI to MNI152 1mm space and downsample to MNI152 2mm space
#########################################################################
foreach runfolder ($bold)
	pushd $runfolder
	echo "Run: $runfolder" |& tee -a $LF
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	echo $BOLD
	
	# if no smooth operation (sm=0), the final output is ${subject}_bld${runfolder}${BOLD_stem}_FS1mm_MNI1mm_MNI2mm.nii.gz
	# if smooth (sm!=0), the final output is ${subject}_bld${runfolder}${BOLD_stem}_FS1mm_MNI1mm_MNI2mm_sm${sm}.nii.gz
	set output_MNI2mm = $volfolder/${BOLD}_FS1mm_MNI1mm_MNI2mm.nii.gz
	if($sm > 0) then
		set final_output = $volfolder/${BOLD}_FS1mm_MNI1mm_MNI2mm_sm${sm}.nii.gz
	else
		set final_output = $output_MNI2mm
	endif

	# if final mask is applied, the final output should have _finalmask postfix
	if( $?final_mask ) then
		set final_output_mask = `basename $final_output .nii.gz`
		set final_output_mask = "$volfolder/${final_output_mask}_finalmask.nii.gz"
	endif
	
	# check if final output already exists
	if( $?final_mask ) then
		if(-e $final_output_mask) then
			echo "[native2mni]: final output: $final_output_mask already exists." |& tee -a $LF
			popd
			continue
		endif
	else
		if(-e $final_output) then
			echo "[native2mni]: final output: $final_output already exists." |& tee -a $LF
			popd
			continue
		endif
	endif
	
	# check if registraton file exists
	set regfile = ${subject}_bld${runfolder}${reg_stem}.dat
	if(! -e $regfile) then
		echo "ERROR: registration file $regfile not exists." |& tee -a $LF
		exit 1;
	endif
	
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
	
	###################################
	### project to MNI152 1mm space
	###################################
	echo "======== Project $runfolder to MNI152 1mm space ========" |& tee -a $LF
	set fcount = 0;
	while($fcount < $nframes)
		set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`
		
		set input = $frame_dir/orig_frames${fcount_str}.nii.gz
		set output = $frame_dir/${fcount_str}_MNI1mm.nii.gz
		if(-e $output) then
			echo "    [native2mni]: $output already exists." |& tee -a $LF
		else
			set cmd = (CBIG_vol2vol_m3z.csh -src-id $anat_s -src-dir $anat_dir -targ-id $MNI_ref_id -targ-dir $MNI_ref_dir -in $input -out $output -reg $regfile -no-cleanup)
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
			if(-e $output) then
				echo "    [native2mni]: projection to $output finished." |& tee -a $LF
			else
				echo "    ERROR: projection to $output failed." |& tee -a $LF
				exit 1;
			endif
		endif
		
		@ fcount = $fcount + 1
	end
	echo "======== Projection of $runfolder to MNI152 1mm space finished. ========" |& tee -a $LF
	echo "" |& tee -a $LF
	
	###################################
	### downsample to MNI152 2mm space
	###################################
	echo "======== Downsample $runfolder to MNI152 2mm space ========" |& tee -a $LF
	setenv SUBJECTS_DIR $MNI_ref_dir
	set fcount = 0;
	while($fcount < $nframes)
		set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`
		set input = $frame_dir/${fcount_str}_MNI1mm.nii.gz
		set output = $frame_dir/${fcount_str}_MNI1mm_MNI2mm.nii.gz
		
		if(-e $output) then
			echo "[native2mni]: $output already exists." |& tee -a $LF
		else
			set cmd = (mri_vol2vol --mov $input --s $MNI_ref_id --targ $temp_2mm --o $output --regheader --no-save-reg)
			echo $cmd |& tee -a $LF
			eval $cmd
			if(-e $output) then
				echo "    [native2mni]: downsample to $output finished." |& tee -a $LF
			else
				echo "    ERROR: downsample to $output failed." |& tee -a $LF
			endif
		endif
		
		@ fcount = $fcount + 1
	end
	echo "" |& tee -a $LF
	
	#########################
	### smooth
	#########################
	echo "======== Smooth in MNI152 2mm space with fwhm = $sm ========" |& tee -a $LF
	set fcount = 0;
        
	if($?sm_mask) then
		set inverted_sm_mask = $frame_dir/inverted_sm_mask.nii.gz
        	set cmd = (fslmaths $sm_mask -mul -1 -add 1 -mas $final_mask $inverted_sm_mask)
        	echo $cmd |& tee -a $LF
        	eval $cmd
  	endif

	while($fcount < $nframes)
		set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`
		
		set input = $frame_dir/${fcount_str}_MNI1mm_MNI2mm.nii.gz
		mkdir -p $frame_dir/sm
		set output = $frame_dir/sm/${fcount_str}_MNI1mm_MNI2mm_sm${sm}.nii.gz
		set std = `awk "BEGIN {print ${sm}/2.35482}"`  #Note that fwhm = 2.35482 * std, fslmaths -s is in unit of mm, not voxel.
		if(-e $output) then
			echo "[native2mni]: $output already exists." |& tee -a $LF
		else
			if($?sm_mask) then
			    # if the user passes in a volume sm_mask, the procedure is
			    # 1. smooth volume data within sm_mask
			    # 2. smooth sm_mask within sm_mask
			    # 3. divide smoothed volume by smoothed sm_mask (deal with boundary problem)
                            set tmp1 = $frame_dir/tmp1_${fcount_str}.nii.gz
	  		    set tmp2 = $frame_dir/tmp2_${fcount_str}.nii.gz
	                    set tmp3 = $frame_dir/tmp3_${fcount_str}.nii.gz		
                            set tmp4 = $frame_dir/tmp4_${fcount_str}.nii.gz
                            set input_masksmoothed = $frame_dir/input_masksmoothed_${fcount_str}.nii.gz
                            set input_outsidemask_smoothed = $frame_dir/input_outsidemask_smoothed_${fcount_str}.nii.gz

			    set cmd = (fslmaths $input -mas $sm_mask -s $std -mas $sm_mask $tmp1)
			    echo $cmd |& tee -a $LF
			    eval $cmd
				
			    set cmd = (fslmaths $sm_mask -s $std -mas $sm_mask $tmp2)
		            echo $cmd |& tee -a $LF
			    eval $cmd
				
			    set cmd = (fslmaths $tmp1 -div $tmp2 $input_masksmoothed)
		            echo $cmd |& tee -a $LF
			    eval $cmd

                            #4. smooth outside of the mask, but inside the brain			    
			    set cmd = (fslmaths $input -mas $inverted_sm_mask -s $std -mas $inverted_sm_mask $tmp3)
			    echo $cmd |& tee -a $LF
			    eval $cmd

			    set cmd = (fslmaths $inverted_sm_mask -s $std -mas $inverted_sm_mask $tmp4)
			    echo $cmd |& tee -a $LF
			    eval $cmd

			    set cmd = (fslmaths $tmp3 -div $tmp4 $input_outsidemask_smoothed)
			    echo $cmd |& tee -a $LF
			    eval $cmd

                            #5. combine externally smoothed data with rest of the data
			    set cmd = (fslmaths $input_masksmoothed -add $input_outsidemask_smoothed $output)
                            echo $cmd |& tee -a $LF
	                    eval $cmd
				
			else
			    set cmd = (fslmaths $input -s $std $output)
			    echo $cmd |& tee -a $LF
			    eval $cmd
			endif
			
			if(-e $output) then
			    echo "[native2mni]: smooth to $output finished." |& tee -a $LF
			else
			    echo "ERROR: smooth to $output failed." |& tee -a $LF
			    exit 1;
			endif
		endif
		
		@ fcount = $fcount + 1
	end
	echo "" |& tee -a $LF
	
	######################
	### combine frames
	######################
	echo "======== Combine frames for $runfolder ========" |& tee -a $LF
	if($sm > 0) then
		set cmd = (fslmerge -t $final_output $frame_dir/sm/*_MNI1mm_MNI2mm_sm${sm}.nii.gz)
	else
		set cmd = (fslmerge -t $final_output $frame_dir/*_MNI1mm_MNI2mm.nii.gz)
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
	
	#########################
	### apply final mask
	#########################
	if( $?final_mask ) then
		echo "======== Applying final mask for $runfolder ========" |& tee -a $LF
		set input = $final_output
		set output = $final_output_mask
		
		if( -e $output ) then
			echo "[native2mni]: $output already exists" |& tee -a $LF
		else		
			set cmd = (fslmaths $input -mas ${final_mask} $output)
			echo $cmd
			eval $cmd

			if( $nocleanup == 0 ) then
				rm ${input}
			endif
			
			if( -e $output ) then
				echo "[native2mni]: Applying final mask finished. The output is $output" |& tee -a $LF
			else
				echo "ERROR: Applying final mask failed." |& tee -a $LF
			endif			
		endif
	endif
	
	
	#########################
	### clean up
	#########################
	if( $nocleanup == 0 ) then
		rm -R $frame_dir
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
	pushd ${CBIG_CODE_DIR}
        git log -1 -- ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2mni.csh >> $LF
        popd
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
		#subject name (required)
		case "-s":
			if ( $#argv == 0 ) goto arg1err;
			set subject = $argv[1]; shift;
			breaksw	
		
		#path to subject's folder (required)
		case "-d":
			if ( $#argv == 0 ) goto arg1err;
			set sub_dir = $argv[1]; shift;
			breaksw
			
		#anatomical directory (required)
		case "-anat_d":
			if ( $#argv == 0 ) goto arg1err;
			set anat_dir = $argv[1]; shift;
			breaksw
		
		#anatomical data (required)
		case "-anat_s":
			if ( $#argv == 0 ) goto arg1err;
			set anat_s = $argv[1]; shift;
			breaksw
		
		# bold number, like "002 003"	
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set bold = ($argv[1]); shift;
			breaksw
			
		#MNI space smooth param
		case "-sm":
			if ( $#argv == 0 ) goto arg1err;
			set sm = $argv[1]; shift;
			breaksw
			
		#smooth mask
		case "-sm_mask":
			if ( $#argv == 0 ) goto arg1err;
			set sm_mask = $argv[1]; shift;
			breaksw
		
		#BOLDbase_suffix 
		case "-BOLD_stem":
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = "$argv[1]"; shift;
			breaksw
			
		#reg suffix (default rest_skip*_stc_mc_reg)
		case "-REG_stem":
			if ( $#argv == 0 ) goto arg1err;
			set reg_stem = "$argv[1]"; shift;
			breaksw
		
		# downsample flag, determine which MNI 2mm space the user wants to use
		# "FSL_MNI_FS_2mm" for FreeSurfer 2mm space (128*128*128); "FSL_MNI_2mm" for FSL MNI 2mm space (91*109*91)
		case "-down":
			if ( $#argv == 0 ) goto arg1err;
			set down = $argv[1]; shift;
			breaksw
			
		# final mask (a loose MNI brain mask)
		case "-final_mask":
			if ( $#argv == 0 ) goto arg1err;
			set final_mask = $argv[1]; shift;
			breaksw
			
		# if nocleanup is turned on, it won't delete the folders with single frames, 
                # otherwise these folders will be removed
		case "-nocleanup":
			set nocleanup = 1
			breaksw
		
		#update results, if exist then do not generate
		case "-force":
			set force = 1;
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

# ${down} determine which downsample template will be used
# FSL_MNI_FS_2mm (128x128x128) or FSL_MNI_2mm (91*109*91)
if( $?down ) then
	if( $down == "FSL_MNI_FS_2mm" ) then
		set temp_2mm = ${FS_temp_2mm}
	else if ( $down == "FSL_MNI_2mm" ) then
		set temp_2mm = ${MNI_temp_2mm}
	else
		echo "ERROR: Wrong input for -down. Please choose from FSL_MNI_FS_2mm and FSL_MNI_2mm. Current down = $down"
		exit 1;
	endif
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
	CBIG_preproc_native2mni.csh
	
DESCRIPTION:
	Project data from fMRI native space to MNI152 1mm space, downsample to MNI152 2mm space, and smooth.
	
	This function does the following steps:
	    1. Project anatomical data to MNI152 1mm space for check purpose
	    2. Project fMRI volumes to MNI152 1mm space
	    3. Downsample the projected data to MNI152 2mm space (either FS 2mm or FSL 2mm, specified by -down)
	    4. Smooth the downsampled data. The user can specify a mask for smoothing.
	    5. Apply the final mask, if any.

	Projection step uses mri_vol2vol to first project volume data to FreeSurfer nonlinear space and 
	then project it to MNI152 1mm space. For more details, please refer to the command 
	'CBIG_vol2vol_m3z.csh'.

	For smoothing step, the users can pass in a mask (using -sm_mask flag). If the users specify a mask, 
	then the program will smooth the volume within the mask, smooth the mask, and divide the smoothed 
	volume by the smoothed mask (to deal with boundary problems).

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
	-sm         sm        : smooth fwhm (mm), e.g. 6. Even if the user does not use this option, 
	                        the scripts will automatically smooth the volume by 6mm. If the user 
	                        does not want to do smoothing, he/she needs to pass in -sm 0.
	-sm_mask    sm_mask   : mask for smoothing (e.g. a grey matter mask in MNI152 2mm). An example 
	                        of the smooth mask is: 
	                        ${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152_masks/SubcorticalLooseMask_MNI1mm_sm6_MNI2mm_bin0.2.nii.gz
	                        If <sm_mask> is not passed in, the smoothing step will smooth everything
	                        by the FWHM as specified by -sm flag.
	-down       down      : downsample space, choose from FSL_MNI_FS_2mm (size: 128 x 128 x 128)
	                        or FSL_MNI_2mm (size: 91 x 109 x 91). Default is FSL_MNI_2mm
	-final_mask final_mask: a loose mask in MNI152 2mm space applied to the volume in the
	                        last step to reduce space. An example of the final mask is:
	                        ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
	-nocleanup            : do not remove intermediate results.
	
OUTPUTS:
	(1) The anatomical volume projected to MNI152 1mm space (for checking):
	    <sub_dir>/<subject>/vol/norm_MNI152_1mm.nii.gz
	
	(2) The fMRI volume after projection (in MNI152 1mm space):
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI1mm.nii.gz
	
	(3) The fMRI volume after downsample (in MNI152 2mm space):
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI1mm_MNI2mm.nii.gz
	
	(4) The fMRI volume after smooth:
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI1mm_MNI2mm_sm<sm>.nii.gz
	
	(5) The fMRI volume after applying the final mask (if -final_mask is used):
	    <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem>_MNI1mm_MNI2mm_sm<sm>_finalmask.nii.gz
	
	(6) Intermediate results (if -nocleanup is passed in):
	    Folders <sub_dir>/<subject>/vol/<subject>_bld<run_number><BOLD_stem> including all intermediate 
	    fMRI single frames. (The projection, downsample, and smoothing are done frame by frame to save 
	    memory)
	    <sub_dir>/<subject>/vol/FStmp.norm_MNI152_1mm.nii.gz is the anatomical volume project to FreeSurfer 
	    nonlinear space.

Example:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2mni.csh 
	-s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS -anat_d ~/storage/sMRI_preprocess 
	-bld '002 003' -BOLD_stem _rest_stc_mc_cen_resid_lp0.08 -REG_stem _rest_stc_mc_reg -down FSL_MNI_2mm -sm 6 -sm_mask
        ${CBIG_CODE_DIR}/data/templates/volume/FSL_MNI152_masks/SubcorticalLooseMask_MNI1mm_sm6_MNI2mm_bin0.2.nii.gz
	-final_mask ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
