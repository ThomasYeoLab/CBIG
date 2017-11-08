#!/bin/csh -f

#############################################
# Detrend the data and do bandpass filtering
#############################################
# In this script, we: 
# 1) demean and remove the linear trend (optional)
# 2) use GLM regression to do bandpass filtering, check function CBIG_bpss_by_regression.m
# 3) recover the linear trend of the original input (optional) 
# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_bandpass_fft.csh 
#	-s Sub0060_Ses1 -d ~/storage/fMRI_preprocess -bld "002 003"	-BOLD_stem _rest_skip_stc_mc_cen_resid 
#	-low_f 0 -high_f 0.08 -detrend -censor
#############################################
# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#############################################

set VERSION = '$Id: CBIG_preproc_bandpass_fft.csh, v 1.0 2016/06/18 $'

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

set subject = ""
set subject_dir = ""
set zpdbold = ""
set BOLD_stem = ""
set OUTLIER_stem = ""

set detrend = 0 	# Default not demean and remove the linear trend
set retrend = 0 	# Default not recover the linear trend of the original input
set low_f = "" 		# Default not set the low_f cutoff frequency
set high_f = "" 	# Default not set the high_f cutoff frequency
set censor = 0		# Default use all frames to do GLM regression
set force = 0		# Default if file exist, skip the step

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

cd $subject_dir/$subject

#############################################
# BOLD Information
#############################################

if (! -e logs) then
    mkdir -p logs
endif
set LF = $subject_dir/$subject/logs/CBIG_preproc_bandpass.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[Bandpass]: logfile = $LF"
echo "bandpass" >> $LF
echo "[CMD]: CBIG_preproc_bandpass.csh $cmdline"   >>$LF

# boldfolder: directory of /bold
set boldfolder = "$subject_dir/$subject/bold"
echo "[Bandpass]: boldfolder = $boldfolder" |& tee -a $LF
echo "[Bandpass]: zpdbold = $zpdbold" |& tee -a $LF 

# go to boldfolder
cd $boldfolder

# check if matlab exists
set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
if ($status) then
    echo "ERROR: could not find matlab"
    exit 1;
endif


#############################################
# Use matlab to do bandpass for each run
#############################################

echo "=======================Bandpass each run =======================" |& tee -a $LF
foreach curr_bold ($zpdbold)

	# Get censor file
	if ( $censor == 1 ) then
		set censor_file = "$subject_dir/$subject/qc/$subject"_bld"${curr_bold}${OUTLIER_stem}"
	else 
		set censor_file = ""
	endif

	# get fMRI file
	pushd $curr_bold
	set boldfile = $subject"_bld"$curr_bold$BOLD_stem
	echo "[Bandpass]: boldfile = $boldfile" |& tee -a $LF

	# bandpass
	if ( ( $low_f != "" ) && ( $high_f != "" ) ) then
		if ( (! -e  $boldfile"_bp_"$low_f"_"$high_f".nii.gz") || ($force == 1) ) then
			set fMRI_file = $boldfile".nii.gz"
			set output_file = $boldfile"_bp_"$low_f"_"$high_f".nii.gz"
			$MATLAB -nojvm -nodesktop -nodisplay -nosplash -r "addpath(fullfile('$root_dir','utilities'));CBIG_bandpass_vol '$fMRI_file' '$output_file' '$low_f' '$high_f' '$detrend' '$retrend' '$censor_file';exit " |& tee -a $LF		
		else
			echo "=======================Bandpass has been done!=======================" |& tee -a $LF
		endif
	endif

	popd
end
echo "=======================Bandpass done!=======================" |& tee -a $LF

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_bandpass_fft.csh >> $LF
endif

echo "*********************************************************************" |& tee -a $LF
exit 1;


##########################################
# Parse Arguments 
##########################################	

parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		#subject name
		case "-s":
			if ( $#argv == 0 ) goto arg1err;
			set subject = $argv[1]; shift;
			breaksw	
			
		#path to subject's folder
		case "-d":
			if ( $#argv == 0 ) goto arg1err;
			set subject_dir = $argv[1]; shift;
			breaksw
			
		#zpdbold: all runs under /bold folder
		case "-bld"
			if ( $#argv == 0 ) goto arg1err;
			set zpdbold = ($argv[1]); shift;
			breaksw
			
		#input file stem
		case "-BOLD_stem"
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = $argv[1]; shift;
			breaksw
		
		#stem of file including censor vector 
		case "-OUTLIER_stem"
			if ( $#argv == 0 ) goto arg1err;
			set OUTLIER_stem = $argv[1]; shift;
			breaksw

		#set low_f cutoff Hz
		case "-low_f":
			if ( $#argv == 0 ) goto arg1err;
			set low_f = $argv[1]; shift;
			breaksw
			
		#set high_f cutoff Hz
		case "-high_f":
			if ( $#argv == 0 ) goto arg1err;
			set high_f = $argv[1]; shift;
			breaksw
			
		#demean and remove the linear trend (optional)
		case "-detrend":
			set detrend = 1;
			breaksw	
			
		#recover the linear trend of the original input (optional) 
		case "-retrend":
			set retrend = 1;
			breaksw	
		
		#remove censored frames and do glm regression 
		case "-censor":
			set censor = 1;
			breaksw		
		
		#update results, if exist then overwrite
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


##########################################
# Check Parameters
##########################################

check_params:
if ( "$subject" == "" ) then
	echo "ERROR: subject not specified"
	exit 1;
endif
if ( "$subject_dir" == "" ) then
	echo "ERROR: path to subject folder not specified"
	exit 1;
endif
if ( "$zpdbold" == "" ) then
	echo "ERROR: bold run not specified"
	exit 1;
endif
if ( "$BOLD_stem" == "" ) then
	echo "ERROR: input file stem not specified"
	exit 1;
endif
if ( "$OUTLIER_stem" == "" ) then
	echo "ERROR: input file stem not specified"
	exit 1;
endif

if ( $low_f == "" )  then
	echo "ERROR: low_f low frequency not specified"	
	exit 1;
endif		

if ( $high_f == "" )  then
	echo "ERROR: high_f high frequency not specified"	
	exit 1;
endif		

if ( ( $low_f != "" ) && ( $high_f != "" ) ) then
	if ( `(echo "if($low_f > $high_f) 1" | bc)`) then
        echo "ERROR: low_frequency $low_f > high_frequency $high_f, bandstop is not allowed"
        exit 1;
	endif
endif

goto check_params_return;

##########################################
# ERROR message
##########################################
		
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1
  
#####################################
# Help
#####################################
BEGINHELP

NAME:
	CBIG_preproc_bandpass_fft.csh

DESCRIPTION:
	Bandpass filtering by regression method.
	
	This function does
	1) demean and remove the linear trend (optional) 
	2) use regression to do bandpass filtering, [low_f high_f] is included the passband window 
	3) recover the linear trend of the original input (optional) 

REQUIRED ARGUMENTS:
	-s  <subject_id>             : subject's ID
	-d  <subject_dir>            : absolute path to <subject_id>. All preprocessed data of this subject 
	                               are assumed to be at <subject_dir>/<subject_id>.
	-bld  <bold_runs>            : all bold run numbers of this subject. Each number must be three 
	                               digits. If this subject has multiple runs, please use space as
	                               delimiter between two runs (e.g. -bld "002 003"). NOTE: quote sign 
	                               is necessary.
	-BOLD_stem  <BOLD_stem>      : the stem of input file. E.g. if input file is 
	                               Sub0001_bld002_rest_skip_stc_mc_resid_cen_FDRMS0.2_DVARS50.nii.gz, 
	                               then <BOLD_stem> = _rest_skip_stc_mc_resid_cen_FDRMS0.2_DVARS50. This 
	                               file is assumed to be stored in <subject_dir>/<subject_id>/bold/<run_number>.
	-OUTLIER_stem <OUTLIER_stem> : the stem of outlier text file. The number of lines of this file is 
	                               equal to the total number of frames. Each line is a number of 0 or 1 
	                               (0-censored, 1-keep). E.g. if the outlier file is 
	                               Sub0001_bld002_FDRMS0.2_DVARS50_motion_outliers.txt, then <OUTLIER_stem> 
	                               = _FDRMS0.2_DVARS50_motion_outliers.txt. This file is assumed to be 
	                               stored in <subject_dir>/<subject>/qc.
	-low_f <low_f>               : low frequency, the boundary low_f will be included in passband window
	-high_f <high_f>             : high frequency, the boundary high_f will be included in passband window

OPTIONAL ARGUMENTS:
	-detrend                     : demean and remove the linear trend before bandpass filtering (optional)
	-retrend                     : add back the linear trend after bandpass filtering (optional), retrend 
	                               is ignored if detrend is off
	-censor                      : fit beta coefficients in regression without the censored frames, then 
	                               apply the beta to all frames.
	-force                       : update results, if exist then overwrite
	-help                        : help
	-version                     : version
	
OUTPUTS:
	The NIFTI volume after bandpass filtering:
	<sub_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_bp_<low_f>_<high_f>.nii.gz

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_bandpass_fft.csh 
	-s Sub0060_Ses1 -d ~/storage/fMRI_preprocess -bld "002 003"	-BOLD_stem _rest_skip_stc_mc_cen_resid 
	-low_f 0 -high_f 0.08 -detrend -censor
	
