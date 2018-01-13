#!/bin/csh -f

# This function needs SUBJECTS_DIR to be well set.
#
# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2fsaverage.csh 
#	-s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS -anat_d ~/storage/sMRI_preprocess 
#	-bld '002 003' -proj fsaverage6 -down fsaverage5 -sm 6 -BOLD_stem 
#	_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0_0.08 -REG_stem _rest_skip4_stc_mc_reg
#
# This function does the following three steps:
#     1. Project fMRI volume data to surface (defined by proj_mesh, like fsaverage6) using mri_vol2surf
#     2. Smooth projected data using mri_surf2surf. Since in FreeSurfer 5.3.0, mri_surf2surf with --cortex flag will set vertices within medial wall mask to be zero, we put the medial wall values back after smoothing.
#     3. Downsample smoothed data to downmesh, e.g. fsaverage5, using mri_surf2surf 
#
# Written by Jingwei Li.
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#BOLD: basename of each run input
#boldfolder: directory of /bold
#surffolder: store all surface results (/surf)
#bold: all runs under /bold folder

set subject = ""       # subject ID
set sub_dir = ""       # directory to subjects
set anat_dir = ""      # directory to recon-all folders
set anat_s = ""        # recon-all folders
set bold = ""          # bold numbers, e.g. '002 003'
set BOLD_stem = ""     # BOLD stem, e.g. _rest_stc_mc_cen_FDRMS0.2_DVARS50_resid_lp0.08
set reg_stem = ""      # registration stem, e.g. _rest_stc_mc_reg

set proj_mesh = fsaverage6
set down_mesh = fsaverage5
set sm = 6;

set all_out_stem = ""
set all_out_mesh = ""

########################
# Print help and version
########################
set VERSION = '$Id: CBIG_preproc_native2fsaverage.csh v 1.0 2016/05/29'

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

###############################
# check if matlab exists
###############################
set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
if ($status) then
	echo "ERROR: could not find matlab"
	exit 1;
endif

setenv SUBJECTS_DIR $anat_dir

set currdir = `pwd`
cd $sub_dir/$subject


###################
# create a log file
###################
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_native2fsaverage.log
if( -e ${LF} ) then
	rm $LF
endif
touch $LF
echo "[SURF]: logfile = $LF"
echo "Surface Projection, Downsample & Smooth" >> $LF
echo "[SURF]: $cmdline"   >>$LF


###############################################################################
# specify bold folder and surf folder (surf folder contains all surface results)
###############################################################################
set boldfolder = "$sub_dir/$subject/bold"
set surffolder = "$sub_dir/$subject/surf"
echo "[SURF]: boldfolder = $boldfolder" |& tee -a $LF

pushd $boldfolder

mkdir -p $surffolder

#############################################################################
### fsaverage needs to be in the anat_dir, if not, create links
#############################################################################
if(! -e $anat_dir/$proj_mesh) then
	ln -s $FREESURFER_HOME/subjects/$proj_mesh $anat_dir/$proj_mesh
endif
if(! -e $anat_dir/$down_mesh) then
	ln -s $FREESURFER_HOME/subjects/$down_mesh $anat_dir/$down_mesh
endif

############################################################################
### project data to proj_mesh
############################################################################
echo "============ Project fMRI volumes to surface $proj_mesh ============" |& tee -a $LF
set all_out_stem = "$all_out_stem ${BOLD_stem}_${proj_short}"
set all_out_mesh = ($all_out_mesh $proj_mesh)
foreach runfolder ($bold)
	pushd $runfolder
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	if (! -e $BOLD.nii.gz) then
		echo "ERROR: input file $BOLD.nii.gz not found" |& tee -a $LF
		exit 1;
	endif
	
	# check if registration file exists
	set regfile = ${subject}_bld${runfolder}${reg_stem}.dat
	if(! -e $regfile) then
		echo "ERROR: registration file $regfile not exists." |& tee -a $LF
		exit 1;
	endif
	
	# projection
	foreach hemi (lh rh)
		set output = $surffolder/$hemi.${BOLD}_$proj_short.nii.gz
		if(-e $output) then
			echo "[SURF]: Projection to $hemi.${BOLD}_$proj_short.nii.gz already exist" |& tee -a $LF
		else
			set cmd = (mri_vol2surf --mov ${BOLD}.nii.gz --reg $regfile --hemi  $hemi --projfrac 0.5 --trgsubject $proj_mesh --o $output --reshape --interp trilinear)
			echo $cmd |& tee -a $LF
			eval $cmd
			
			if(-e $output) then
				echo "[SURF]: Projection to $hemi.${BOLD}_$proj_short.nii.gz finished" |& tee -a $LF
			else
				echo "[SURF]: Projection to $hemi.${BOLD}_$proj_short.nii.gz failed" |& tee -a $LF
			endif
		endif
	end
	
	popd
end
echo "======================= Projection finished. ======================" |& tee -a $LF
echo "" |& tee -a $LF


########################################################
### smooth
########################################################
echo "================== Smooth surface data, fwhm $sm =================" |& tee -a $LF
set all_out_stem = "$all_out_stem ${BOLD_stem}_${proj_short}_sm${sm}"
set all_out_mesh = ($all_out_mesh $proj_mesh)
foreach runfolder ($bold)
	pushd $runfolder
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	foreach hemi (lh rh)
		set input = $surffolder/$hemi.${BOLD}_${proj_short}.nii.gz
		set output = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}.nii.gz
		if(-e $output) then
			echo "[SURF]: Smooth result $output already exists." |& tee -a $LF
		else
			set cmd = (mri_surf2surf --hemi $hemi --s $proj_mesh --sval $input --cortex --fwhm-trg $sm --tval $output --reshape)
			echo $cmd |& tee -a $LF
			eval $cmd
			
			### mri_surf2surf with --cortex flag will make medial wall vertices to be zeros
			### add removed medial wall values back
			echo "===>> Fill in zero medial values " |& tee -a $LF
			set input1 = $surffolder/$hemi.${BOLD}_${proj_short}.nii.gz
			set input2 = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}.nii.gz
			set tmp_output = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}_fillin_medialwall.nii.gz
			
			if( ! -e $input1 || ! -e $input2 ) then
				echo "ERROR: One or two of the input data do not exist!"
				echo "Input1: $input1"
				echo "Input2: $input2"
				exit 1;
			endif
			
			$MATLAB -nodesktop -nodisplay -nosplash -r "addpath(fullfile('$root_dir', 'utilities'));CBIG_preproc_fsaverage_medialwall_fillin '$hemi' 'fsaverage6' '$input1' '$input2' '$tmp_output';exit" |& tee -a $LF
			
			if ( ! -e $tmp_output ) then
				echo "ERROR: Fill in medial wall failed. ${tmp_output} is not produced." |& tee -a $LF
				exit 1
			else
				echo "[SURF]: ${tmp_output} is successfully produced." |& tee -a $LF
			endif

			set cmd = "mv ${tmp_output} ${input2}"
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
			echo "===>> Fill in medial wall values finished " |& tee -a $LF
			
			if(-e $output) then
				echo "[SURF]: Smooth for $output finished." |& tee -a $LF
			else
				echo "[SURF]: Smooth for $output failed." |& tee -a $LF
			endif
		endif
	end
	popd
end
echo "====================== Smoothness finished =====================" |& tee -a $LF
echo "" |& tee -a $LF




###################################################
### downsample
###################################################
echo "==================== Downsample to $down_mesh ==================" |& tee -a $LF
set all_out_stem = "$all_out_stem ${BOLD_stem}_${proj_short}_sm${sm}_${down_short}"
set all_out_mesh = ($all_out_mesh $down_mesh)

if($proj_res < $down_res) then
	echo "ERROR: projection mesh ($proj_mesh) < downsampling mesh ($down_mesh)" |& tee -a $LF
	exit 1;
endif

foreach runfolder ($bold)
	pushd $runfolder
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	foreach hemi (lh rh)
		set input = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}.nii.gz
		set curr_input = $input
		set output = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}_${down_short}.nii.gz
		
		if(-e $output) then
			echo "[SURF]: Downsampling result $output already exists." |& tee -a $LF
		else
			set scale = $proj_res
			if($scale == $down_res) then
				set cmd = (cp $curr_input $output)
				echo $cmd |& tee -a $LF
				eval $cmd
			endif
			
			# downsample scale by scale
			while($scale > $down_res) 
				@ new_scale = $scale - 1
				if($scale == 7) then
					set srcsubject = fsaverage
				else
					set srcsubject = fsaverage$scale
				endif
				
				set trgsubject = fsaverage$new_scale
				
				set cmd = (mri_surf2surf --hemi $hemi --srcsubject $srcsubject --sval $curr_input --nsmooth-in 1 --trgsubject $trgsubject --tval $output --reshape)
				echo $cmd |& tee -a $LF
				eval $cmd
				
				set curr_input = $output
				@ scale = $scale - 1
			end
		endif
	end
	popd
end
echo "=================== Downsampling finished ==================" |& tee -a $LF


#########################
# Set medial values to be NaN
#########################
foreach runfolder ($bold)
	set i = 1;
	foreach out_stem ($all_out_stem)
		set out_mesh = $all_out_mesh[$i];
		# for debugging
		echo "out_stem: $out_stem" |& tee -a $LF
		echo "out_mesh: $out_mesh" |& tee -a $LF
		foreach hemi (lh rh)
			set before_NaN_name = $surffolder/$hemi.${subject}_bld${runfolder}${out_stem}.nii.gz
			set after_NaN_name = $surffolder/$hemi.${subject}_bld${runfolder}${out_stem}_medialwallNaN.nii.gz
			
			$MATLAB -nodesktop -nodisplay -nosplash -r "addpath(fullfile('$root_dir', 'utilities'));CBIG_preproc_set_medialwall_NaN '$hemi' '$out_mesh' '$before_NaN_name' '$after_NaN_name';exit" |& tee -a $LF
			
			if ( -e $after_NaN_name ) then
				echo "[SURF]: $after_NaN_name is successfully generated. Move it to $before_NaN_name." |& tee -a $LF
				mv $after_NaN_name $before_NaN_name
			else
				echo "ERROR: $after_NaN_name is not generated." |& tee -a $LF
				exit 1
			endif
		end
		
		set i = `expr $i + 1`
	end
end


#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2fsaverage.csh >> $LF
endif

echo "***************************************************************" |& tee -a $LF
echo "" |& tee -a $LF


popd

exit 0;

#################################################
##======pass the arguments
#################################################
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
			
		# bold number, e.g. '002 003'
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set bold = ($argv[1]); shift;
			breaksw
			
		#projection mesh (fsaverage, fsaverage6, ...)
		case "-proj"
			if ( $#argv == 0 ) goto arg1err;
			set proj_mesh = $argv[1]; shift;
			breaksw
			
		#downsample mesh (fsaverage, fsaverage6, ...)
		case "-down"
			if ( $#argv == 0 ) goto arg1err;
			set down_mesh = $argv[1]; shift;
			breaksw
			
		#smooth param
		case "-sm"
			if ( $#argv == 0 ) goto arg1err;
			set sm = $argv[1]; shift;
			breaksw
		
		#BOLD stem
		case "-BOLD_stem"
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = "$argv[1]"; shift;
			breaksw
			
		#registration stem
		case "-REG_stem"
			if ( $#argv == 0 ) goto arg2err;
			set reg_stem = "$argv[1]"; shift;
			breaksw
		
		default:
			echo ERROR: Flag $flag unrecognized.
			echo $cmdline
			exit 1
			breaksw
	endsw
end
goto parse_args_return;


################################################################
##======check passed parameters
################################################################
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

# check the format of projection and downsampling mesh
if($proj_mesh != fsaverage & $proj_mesh != fsaverage6 & $proj_mesh != fsaverage5 & $proj_mesh != fsaverage4 ) then
    echo "ERROR: proj_mesh = $proj_mesh is not acceptable (allowable values = fsaverage, fsaverage6, fsaverage5, fsaverage4)"
    exit 1;
endif

if($down_mesh != fsaverage & $down_mesh != fsaverage6 & $down_mesh != fsaverage5 & $down_mesh != fsaverage4 ) then
    echo "ERROR: down_mesh = $down_mesh is not acceptable (allowable values = fsaverage, fsaverage6, fsaverage5, fsaverage4)"
    exit 1;
endif


# projection and downsample resolution
set proj_res = `echo -n $proj_mesh | tail -c -1`
if($proj_res == "e") then
	set proj_res = 7
endif

set down_res = `echo -n $down_mesh | tail -c -1`
if($down_res == "e") then
	set down_res = 7;
endif

set proj_short = "fs$proj_res"
set down_short = "fs$down_res"
				
goto check_params_return;


###########################################		
##======Error message		
###########################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1

arg1err:
  echo "ERROR: flag $flag requires at least one argument"
  exit 1


###########################################
# usage exit
###########################################
BEGINHELP

NAME: 
	CBIG_preproc_native2fsaverage.csh

DESCRIPTION:
	Project data from fMRI native space to surface, smooth, and downsample.
	
	This function does the following three steps:
	    1. Project fMRI volume data to FreeSurfer surface (defined by proj_mesh, like fsaverage6) 
	       using mri_vol2surf.
	    2. Smooth projected data using mri_surf2surf. Since in FreeSurfer 5.3.0, mri_surf2surf with 
	       --cortex flag will set vertices within medial wall mask to be zero, we put the medial wall 
	       values back after smoothing.
	    3. Downsample smoothed data to down_mesh, e.g. fsaverage5, using mri_surf2surf
	    4. Set the medial wall to be NaN for all output files.



REQUIRED ARGUMENTS:
	-s          subject   : fMRI subject id
	-d          sub_dir   : absolute path to <subject>. All preprocessed data of this 
	                        subject are assumed to be stored in <sub_dir>/<subject>
	-anat_s     anat_s    : anatomical subject id (recon-all output subject name)
	-anat_d     anat_dir  : absolute path to <anat_s>. All recon-all outputs of this 
	                        subject are stored in <anat_dir>/<anat_s>.
	-bld        bold      : all bold run numbers of this subject. Each number must be 
	                        three digits. If this subject has multiple runs, use a space 
	                        as delimiter, e.g. '002 003'. NOTE: quote sign is necessary.
	-BOLD_stem  BOLD_stem : stem of input volume, e.g. if the input file name is
	                        Sub0001_Ses1_bld002_rest_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0_0.08.nii.gz,
	                        then <BOLD_stem> = _rest_stc_mc_cen_resid_cen_FDRMS0.2_DVARS50_bp_0_0.08. 
	                        This file is assumed to be stored in <sub_dir>/<subject>/bold/<run_number>.
	-REG_stem   reg_stem  : stem of T1-T2* registration file. E.g. if the registration 
	                        file is Sub0001_Ses1_bld002_rest_skip4_stc_mc_reg.dat, 
	                        then <reg_stem> = _rest_skip4_stc_mc_reg. The registration 
	                        file is assumed to be stored in <sub_dir>/<subject>/bold/<run_number>.

OPTIONAL ARGUMENTS:
	-proj       proj_mesh : projection resolution, e.g. fsaverage6 (default)
	-down       down_mesh : downsample resolution, e.g. fsaverage5 (default)
	-sm         sm        : smooth fwhm (mm), e.g. 6 (default)

OUTPUTS:
	Three NIFTI volumes will be output.
	(1) The volume after projection (e.g. proj_mesh = fsaverage6):
	    <sub_dir>/<subject>/bold/<run_number>/<subject>_bld<run_number><BOLD_stem>_fs6.nii.gz
	(2) The volume after smooth (e.g. proj_mesh = fsaverage6, sm = 6):
	    <sub_dir>/<subject>/bold/<run_number>/<subject>_bld<run_number><BOLD_stem>_fs6_sm6.nii.gz
	(3) The volume after downsample (e.g. proj_mesh = fsaverage6, sm = 6, down_mesh = fsaverage5):
	    <sub_dir>/<subject>/bold/<run_number>/<subject>_bld<run_number><BOLD_stem>_fs6_sm6_fs5.nii.gz

Example:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2fsaverage.csh 
	-s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS -anat_d ~/storage/sMRI_preprocess 
	-bld '002 003' -proj fsaverage6 -down fsaverage5 -sm 6 -BOLD_stem 
	_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0_0.08 -REG_stem _rest_skip4_stc_mc_reg



