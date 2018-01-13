#!/bin/csh -f

# example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_censor.csh 
#	-s Sub001_Ses1 -d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS
#	-anat_d ~/storage/sMRI_preprocess -bld '002 003' -BOLD_stem _rest_stc_mc_resid -REG_stem
#	_rest_stc_mc_reg -OUTLIER_STEM _FDRMS0.2_DVARS50_motion_outliers.txt -nocleanup
#
# This function does the following three steps:
#	1. Create grey matter and whole brain masks (if they do not exist)
#	2. Censor "bad" time points (implementation of Power et al. 2014). Detrend before censoring.
#	3. Censoring QC. Compute the correlation and fractional difference between original signal and final signal. Plot the 0%, 20%, 40%, 60%, 80%, and 100% ranking voxels' time series (correlation and fractional difference respectively)
#
# Reference:
# 1) Power, Jonathan D., et al. "Methods to detect, characterize, and remove motion artifact in resting state fMRI." Neuroimage 84 (2014): 320-341.
#
# Written by Jingwei Li.
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#BOLD: basename of each run input
#boldfolder: directory of /bold
#bold: all runs under /bold folder


set subject = ""           # subject ID
set sub_dir = ""           # directory to subjects
set anat_s = ""            # recon-all folder
set anat_dir = ""          # directory to recon-all folder
set bold = ""              # bold number, like '002 003'
set BOLD_stem = ""         # BOLD stem, e.g. _rest_stc_mc
set reg_stem = ""          # registration stem, e.g. _rest_stc_mc_reg
set outlier_stem = ""      # outlier vector file name
set low_f = ""             # low cut-off frequency
set high_f = ""            # high cut-off frequency

set VERSION = '$Id: CBIG_preproc_censor.csh v 1.0 2016/05/26'

set nocleanup = 0;

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

###############################
# check if matlab exists
###############################
set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
if ($status) then
	echo "ERROR: could not find matlab"
	exit 1;
endif

###############################
# create log file
###############################
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_censor.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[CENSOR]: logfile = $LF"
echo "Censoring / Motion Scrubbing" >> $LF
echo "[CMD]: CBIG_preproc_censor.csh $cmdline"   >>$LF

###############################
# echo bold folder
###############################
set boldfolder = "$sub_dir/$subject/bold"
echo "[CENSOR]: boldfolder = $boldfolder" |& tee -a $LF


###############################
# echo FD and DVARS threshold
###############################
set FD = (`echo $outlier_stem | awk -F "FDRMS" '{print $2}' | awk -F "_" '{print $1}'`)
set DVARS = (`echo $outlier_stem | awk -F "DVARS" '{print $2}' | awk -F "_" '{print $1}'`)

echo "FD = $FD"
echo "DVARS = $DVARS"


##############################
# push to bold folder
##############################
echo "" |& tee -a $LF
pushd $boldfolder


###############################################
# Make whole brain and grey matter mask
###############################################
echo "========================= Make whole brain & Grey Matter Mask =========================" |& tee -a $LF
if( -e $boldfolder/mask/${subject}.func.gm.nii.gz && -e $boldfolder/mask/${subject}.brainmask.bin.nii.gz && \
    -e $boldfolder/mask/${subject}.loosebrainmask.bin.nii.gz ) then
	# if grey matter mask and whole brain mask already exist, do nothing
	echo "[CENSOR]: Grey mask and whole brain mask of subject $subject already exist." |& tee -a $LF
else
	# if one or both of GM and whole brain masks does not exist, create them.
	set cmd = "${root_dir}/CBIG_preproc_create_mask.csh -s $subject -d $sub_dir -anat_s ${anat_s} -anat_d ${anat_dir}"
	set cmd = "$cmd -bld '${bold}' -REG_stem ${reg_stem} -MASK_stem ${BOLD_stem} -whole_brain -gm -loose_whole_brain"
	echo $cmd |& tee -a $LF
	eval $cmd
endif
echo "===================== Making whole brain & Grey Matter Mask finished ===================" |& tee -a $LF
echo "" |& tee -a $LF

###############################################
# Censoring
###############################################
echo "===================== Censoring / Motion Scrubbing ======================" |& tee -a $LF

foreach runfolder ($bold)
	echo "Run: $runfolder" |& tee -a $LF
	pushd $runfolder
	
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	set outlier_file = "../../qc/${subject}"_bld"${runfolder}${outlier_stem}"
	set loose_mask = "$boldfolder/mask/${subject}.loosebrainmask.bin.nii.gz"
	
	if ( "$bandpass_flag" == "0" ) then
		# if do not perform bandpass filtering, the output stems do not contain information of passband
		set output_inter = "${BOLD}_interp_inter_FDRMS${FD}_DVARS${DVARS}.nii.gz"
		set output = ${BOLD}_interp_FDRMS${FD}_DVARS${DVARS}.nii.gz
	else
		# if perform bandpass filtering, the output stems contain information of passband
		set output_inter = "${BOLD}_interp_inter_FDRMS${FD}_DVARS${DVARS}_bp_${low_f}_${high_f}.nii.gz"
		set output = ${BOLD}_interp_FDRMS${FD}_DVARS${DVARS}_bp_${low_f}_${high_f}.nii.gz
	endif
	
	# obtain TR information
	set MRI_info = `mri_info $BOLD.nii.gz`
	set TR = `echo $MRI_info | grep -o 'TR: \(.*\)' | awk -F " " '{print $2}'`
	
	if(-e $output && -e $output_inter) then
		# if intermediate and final outputs exist, do nothing.
		echo "$output and $output_inter already exists." |& tee -a $LF
	else
		# if one or both of the two outputs does not exist, call interpolation function
		if ( "$bandpass_flag" == 0 ) then
			# if do not perform bandpass filtering, do not pass in low_f and high_f
			set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r)
			set cmd = ($cmd '"''addpath(fullfile('"'"${root_dir}"'"',' "'"utilities"'"'))';)
			set cmd = ($cmd 'CBIG_preproc_censor_wrapper('"'"${BOLD}.nii.gz"'"',' "'"${outlier_file}"'"',' "'"$TR"'"',' )
			set cmd = ($cmd "'"${output_inter}"'"',' "'"${output}"'"',' "'"${loose_mask}"'"')'; exit '"')
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
		else
			# if perform bandpass filtering, pass in low_f and high_f
			set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r)
			set cmd = ($cmd '"''addpath(fullfile('"'"${root_dir}"'"',' "'"utilities"'"'))';)
			set cmd = ($cmd 'CBIG_preproc_censor_wrapper('"'"${BOLD}.nii.gz"'"',' "'"${outlier_file}"'"',' "'"$TR"'"',' )
			set cmd = ($cmd "'"${output_inter}"'"',' "'"${output}"'"',' "'"${loose_mask}"'"',' "'"${low_f}"'"',')
			set cmd = ($cmd "'"${high_f}"'"')'; exit '"' )
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
		endif
		
		# check whether censoring interpolation successfully exited.
		if(-e $output) then
			echo "Censoring for $runfolder succeed." |& tee -a $LF
		else
			echo "ERROR: Censoring failed." |& tee -a $LF
			exit 1;
		endif
	endif
	
	popd
end
echo "========================= Censoring finished =======================" |& tee -a $LF
echo "" |& tee -a $LF


#############################################
# Censoring QC
#############################################
echo "========================== Censoring QC ==========================" |& tee -a $LF

foreach runfolder ($bold)
	echo "Run: $runfolder" |& tee -a $LF
	pushd $runfolder
	
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	set outlier_file = "../../qc/${subject}"_bld"${runfolder}${outlier_stem}"
	
	if ( "$bandpass_flag" == "0" ) then
		# if do not perform bandpass filtering, the output stems do not contain passband information
		set censor_inter = "${BOLD}_interp_inter_FDRMS${FD}_DVARS${DVARS}.nii.gz"
		set censor_final = ${BOLD}_interp_FDRMS${FD}_DVARS${DVARS}.nii.gz
	else
		# if perform bandpass filtering, the output stems contain passband information
		set censor_inter = "${BOLD}_interp_inter_FDRMS${FD}_DVARS${DVARS}_bp_${low_f}_${high_f}.nii.gz"
		set censor_final = ${BOLD}_interp_FDRMS${FD}_DVARS${DVARS}_bp_${low_f}_${high_f}.nii.gz
	endif
	
	# perform censoring QC
	set qc_out_dir = "${sub_dir}/${subject}/qc/censor_interp"
	mkdir -p $qc_out_dir
	if ( -e ${censor_inter} && -e ${censor_final} ) then
		$MATLAB -nodesktop -nodisplay -nosplash -r "addpath(fullfile('${root_dir}', 'utilities')); CBIG_preproc_CensorQC('${qc_out_dir}', '$subject', '$runfolder', '${BOLD}.nii.gz', '${censor_inter}', '${censor_final}', '${sub_dir}/${subject}/bold/mask/${subject}.brainmask.bin.nii.gz', '${sub_dir}/${subject}/bold/mask/${subject}.func.gm.nii.gz', '${outlier_file}'); exit" |& tee -a $LF
	else
		echo "ERROR: One or two of the intermediate censoring result and final censoring result do not exist. Censoring QC cannot execute." |& tee -a $LF
	endif
	
	# if -nocleanup is not turned on, delete the intermediate volume.
	if( $nocleanup == 0 ) then
		rm $censor_inter
	endif
	
	popd
end
echo "========================= Censoring QC finished =======================" |& tee -a $LF

popd

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_censor.csh >> $LF
endif

exit 0;


##################################
##======pass the arguments
##################################
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
			set sub_dir = $argv[1]; shift;
			breaksw

		case "-anat_s":
			if ( $#argv == 0 ) goto arg1err;
			set anat_s = $argv[1]; shift;
			breaksw

		case "-anat_d":
			if ( $#argv == 0 ) goto arg1err;
			set anat_dir = $argv[1]; shift;
			breaksw
                        
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set bold = ($argv[1]); shift;
			breaksw
			
		#BOLDbase_suffix
		case "-BOLD_stem":
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = "$argv[1]"; shift;
			breaksw
			
		case "-REG_stem":
			if ( $#argv == 0 ) goto arg1err;
			set reg_stem = "$argv[1]"; shift;
			breaksw
			
		case "-OUTLIER_stem":
			if ( $#argv == 0 ) goto arg1err;
			set outlier_stem = $argv[1]; shift;
			breaksw
			
		case "-low_f":
			if ( $#argv == 0 ) goto arg1err;
			set low_f = $argv[1]; shift;
			breaksw
			
		case "-high_f":
			if ( $#argv == 0 ) goto arg1err;
			set high_f = $argv[1]; shift;
			breaksw
			
		case "-nocleanup":
			set nocleanup = 1;
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
if ( "$subject" == "" ) then
	echo "ERROR: subject not specified"
	exit 1;
endif
 
if ( "$sub_dir" == "" ) then
	echo "ERROR: path to subject folder not specified"
	exit 1;
endif

if ( "$anat_s" == "" ) then
	echo "ERROR: anatomical id not specified"
	exit 1;
endif

if ( "$anat_dir" == "" ) then
	echo "ERROR: path to anatomical data not specified"
	exit 1;
endif

if ( "$bold" == "" ) then
        echo "ERROR: bold number not specified"
        exit 1;
endif

if ( "$BOLD_stem" == "" ) then
	echo "ERROR: BOLD stem not specified"
	exit 1;
endif

if ( "$reg_stem" == "" ) then
	echo "ERROR: reg stem not specified. For censoring step, we need to use registration information to create gray matter mask. Please set CBIG_preproc_bbregister before CBIG_preproc_censor in your configuration file."
	exit 1;
endif

if ( "$outlier_stem" == "" ) then
	echo "ERROR: outlier stem not specified"
	exit 1;
endif

if ( "$low_f" == "" && "$high_f" == "" ) then
	set bandpass_flag = 0;
else if ( "$low_f" != "" && "$high_f" != "" ) then
	set bandpass_flag = 1;
else
	echo "low_f or high_f does not exist! Please check the input arguments."
	exit 1;
endif

			
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
	CBIG_preproc_censor.csh

DESCRIPTION:
	Censoring interpolation. Bandpass filtering can be optionally integrated in the 
	interpolation, but it is not recommended. To execute this function, we assume 
	that the T1-T2* registration of this subject has been finished.
	
	This function does the following three steps:
		1. Create grey matter mask and whole brain mask (if the masks do not exist)
		2. Censoring and interpolation of 'bad' time points (implementation of Power et al. 2014).  
		   Linear trend ax+b of timeseries will be removed before censoring. The trend is computed 
		   from uncensored (low motion) frames. (Bandpass filtering can be optionally integrated, 
		   but it is not recommended.)
		3. Censoring QC. Compute the correlation and absolute fractional difference between original 
		   signal and final signal. Save out these two measurements in nifti files and corresponding 
		   statistics (max, min, mean, median) in text files. Plot the 0%, 20%, 40%, 60%, 80%, and 
		   100% ranking voxels' time series (correlation and absolute fractional difference respectively).

	If low_f and high_f are both passed in, this function will perform bandpass filtering. And the 
	final result is the interpolation result within passband. But this usage of integrated interpolation 
	and bandpass filtering is not recommended. We recommend the users to follow the sequence in Power et 
	al. 2014 (interpolate censored frames by this function, then apply bandpass filtering by other 
	functions).
	If the two frequencies are not passed in, this function will only do interpolation and then 
	replace the uncensored ("good") frames with the original signal. If -nocleanup option is used, 
	no matter low_f and high_f are passed in or not, this function will output an intermediate volume 
	which is the interpolation result without bandpass filtering and without replacement.

REQUIRED ARGUMENTS:
	-s             subject      : fMRI subject id
	-d             sub_dir      : absolute path to <subject>. All preprocessed data of this 
	                              subject are assumed to be stored in <sub_dir>/<subject>.
	-anat_s        anat_s       : subject's anatomical id (recon-all output subject name)
	-anat_d        anat_dir     : absolute path to <anat_s>. All recon-all outputs of this 
	                              subject are stored in <anat_dir>/<anat_s>.
	-bld           bold         : all bold run numbers of this subject. Each number must has 
	                              three digits. If this subject has multiple runs, use a space 
	                              as delimiter, e.g. '002 003'. NOTE: quote sign is necessary.
	-BOLD_stem     BOLD_stem    : stem of input volume. E.g. if the input file name is
	                              Sub0001_Ses1_bld002_rest_skip4_stc_mc_resid.nii.gz, then 
	                              <BOLD_stem> = _rest_stc_mc_resid. This input file is assumed to be 
	                              stored in <sub_dir>/<subject>/bold/<run_number>.
	-REG_stem      reg_stem     : stem of T1-T2* registration file. E.g. if the registration 
	                              file is Sub0001_Ses1_bld002_rest_stc_mc_reg.dat, then 
	                              <reg_stem> = _rest_stc_mc_reg. The registration file is 
	                              assumed to be stored in <sub_dir>/<subject>/bold/<run_number>.
	-OUTLIER_stem  outlier_stem : outlier text file stem. The number of lines in this file is 
	                              equal to the total number of frames. Each line is a number of 
	                              0 or 1 (0-censored, 1-keep). E.g. if the outlier file is 
	                              Sub0001_Ses1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt, 
	                              then <outlier_stem> = _FDRMS0.2_DVARS50_motion_outliers.txt. 
	                              The outlier file is assumed to be stored in <sub_dir>/<subject>/qc.
OPTIONAL ARGUMENTS:
	-low_f         low_f        : the low cut-off frequency if bandpass filtering is performed.
	                              E.g., if the passband is [0.009, 0.08], then low_f is 0.009.
	                              Note: the cut-off frequency is also included in the passband. 
	                              However, since this function always to demean on uncensored 
	                              frames and never add it back, if low_f = 0, it is not included 
	                              in passband.
	-high_f        high_f       : the high cut-off frequency if bandpass filtering is performed.
	                              E.g., if the passband is [0.009, 0.08], then high_f is 0.08.
	                              Note: the cut-off frequency is also included in the passband.
	-nocleanup                  : do not remove intermediate result (interpolated volume without 
	                              bandpass filtering and without replacement.

OUTPUTS:
	1. If -nocleanup flag is passed in, this function will output two BOLD volumes:
	   (a) The final output volume (either with interpolation + replacement, or with bandpass + interpolation)
	   <sub_dir>/<subject>/bold/<run_number>/<subject>_bld<run_number><BOLD_stem>_interp_FD<FD_thres>_DVARS<DVARS_thres>.nii.gz
	   
	   (b) The intermediate output volume (all frames are interpolated from input uncensored frames)
	   <sub_dir>/<subject>/bold/<run_number>/<subject>_bld<run_number><BOLD_stem>_interp_inter_FD<FD_thres>_DVARS<DVARS_thres>.nii.gz
	   
	   If -nocleanup is not passed in, only (a) will be output.
	   
	2. Censoring QC plots
	   The voxels are sorted by correlation or absolute fractional difference (both in ascending 
	   ordering) between the original signal and the final censored signal. The voxels with ranking 
	   0%, 20%, 40%, 60%, 80%, and 100% are picked, and three signals of each of these voxels are 
	   plotted, which are: the original signal, the recovered signal (intermediate result), and the 
	   final signal. The procedure are repeated both within whole brain mask and grey matter mask. 
	   The following are the figures that users can check: 
	   
	   The three signals of 6 voxels ranked by correlation within grey matter mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_corr_GM_6plots.jpg
	   
	   The three signals of 6 voxels ranked by absolute fractional difference within grey matter mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_FracDiff_GM_6plots.jpg
	   
	   The three signals of 6 voxels ranked by correlation within whole brain mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_corr_whole_6plots.jpg
	   
	   The three signals of 6 voxels ranked by absolute fractional difference within whole brain mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_FracDiff_whole_6plots.jpg
	
	3. Censoring QC volumes 
	   Besides the censoring QC plots, the correlation and absolute fractional difference of all voxels 
	   will be saved out. If the users want to further know the most different regions before and after 
	   censoring interpolation, they can visualize:
	   
	   Correlation within grey matter mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_corr_GM.nii.gz
	   
	   Absolute fractional difference within grey matter mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_FracDiff_GM.nii.gz
	   
	   Correlation within whole brain mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_corr_whole.nii.gz
	   
	   Absolute fractional difference within whole brain mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_FracDiff_GM.nii.gz
	   
	4. Censoring QC text files
	   Some statistics of the correlation and absolute fractional difference are output. The users can 
	   check the min, max, mean, and median of the correlation and absolute fractional difference. The 
	   text files are:
	   
	   Min, max, mean, and median correlation within grey matter mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_corr_GM.txt
	   
	   Min, max, mean, and median absolute fractional difference within grey matter mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_FracDiff_GM.txt
	   
	   Min, max, mean, and median correlation within whole brain mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_corr_whole.txt
	   
	   Min, max, mean, and median absolute fractional difference within whole brain mask:
	   <sub_dir>/<subject>/qc/<subject>_bld<run_number>_interp_FracDiff_GM.txt

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_censor.csh 
	-s Sub001_Ses1 -d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS
	-anat_d ~/storage/sMRI_preprocess -bld '002 003' -BOLD_stem _rest_stc_mc_resid -REG_stem
	_rest_stc_mc_reg -OUTLIER_STEM _FDRMS0.2_DVARS50_motion_outliers.txt -nocleanup


REFERENCES:
	1) Power, Jonathan D., et al. "Methods to detect, characterize, and remove motion artifact in 
	resting state fMRI." Neuroimage 84 (2014): 320-341.


