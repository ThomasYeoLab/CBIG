#! /bin/csh -f

#############################################
# Skip first several frames
#############################################
# AUTHOR ####################################
# RU(BY) KONG 
# 2016/06/09 
# LYU XINGYU
# 2022/01/14
#############################################
#############################################
# In this script, we: 
# 1) skip first several frames of fMRI data; In the case of multi-echo acquisition, 
#    we will skip the same number of first few frames of each echo

# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_skip.csh 
#	-s Sub0033_Ses1 -d ~/storage/fMRI_preprocess -bld '002 003' -BOLD_stem _rest -skip 4
#############################################
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_skip.csh, v 1.0 2016/06/18 $'

set n = `echo $argv | grep -e -help | wc -l`

# if there is -help option 
if( $n != 0 ) then
	echo $VERSION
	# print help	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
	exit 0;
endif

# if there is no arguments
if( $#argv == 0 ) then
	echo $VERSION
	# print help	
	cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
	echo "WARNING: No input arguments. See above for a list of available input arguments."
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
set skip = 4 #Default skip first 4 frames
set echo_number = 1	#echo number default to be 1
set echo_stem = ""
set force = 0 #Default if file exist, skip the step

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

cd $subject_dir/$subject

#############################################
# BOLD Information
#############################################

if (! -e logs) then
     mkdir -p logs
endif
set LF = $subject_dir/$subject/logs/CBIG_preproc_skip.log
if( -e ${LF} ) then
	rm $LF
endif
touch $LF
echo "[SKIP]: logfile = $LF"
echo "Skip Frames" >> $LF
echo "[CMD]: CBIG_preproc_skip.csh $cmdline"   >>$LF

set boldfolder = "$subject_dir/$subject/bold"
echo "[SKIP]: boldfolder = $boldfolder" |& tee -a $LF
echo "[SKIP]: zpdbold = $zpdbold" |& tee -a $LF

cd $boldfolder

#############################################
# Skip first several frames
#############################################

echo "=======================Skip frames=======================" |& tee -a $LF
foreach curr_bold ($zpdbold)
	pushd $curr_bold
	@ j=1
	while ( $j <= $echo_number)
		if ( $echo_number > 1 ) then 
			set echo_stem = _e$j
		endif
		set boldfile = $subject"_bld${curr_bold}${echo_stem}$BOLD_stem"
		if ( (! -e $boldfile"_skip$skip.nii.gz") || ($force == 1) ) then
			echo "[SKIP]: boldfile = $boldfile" |& tee -a $LF
			@ numof_tps = `fslnvols $boldfile` - $skip
			echo "[SKIP]: Deleting first $skip frames (fslroi) from $boldfile" |& tee -a $LF
			fslroi $boldfile $boldfile"_skip$skip" $skip $numof_tps |& tee -a $LF
			echo "[SKIP]: There are $numof_tps frames after skip $skip frames, output is \
			 $boldfile'_skip$skip.nii.gz'" |& tee -a $LF
			# change file name to be consistant if echo number is 1. 
			if ($echo_number == 1) then
				set cmd = "rsync $boldfile'_skip$skip.nii.gz' $subject"_bld$curr_bold$BOLD_stem"_skip$skip.nii.gz"
				eval $cmd
			endif
		else
			echo "[SKIP]: $boldfile'_skip$skip.nii.gz' already exists!"
		endif

		@ j++
	end
	popd
end
echo "=======================Skip done!=======================" |& tee -a $LF
echo "" |& tee -a $LF

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	pushd ${CBIG_CODE_DIR}
	git log -1 -- ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_skip.csh >> $LF
	popd
endif

echo "****************************************************************" |& tee -a $LF
exit 0;

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
			
		#bold run number	
		case "-bld":
			if ($#argv == 0) goto arg1err;
			set zpdbold = ($argv[1]); shift;
			breaksw
		
		#input file stem	
		case "-BOLD_stem":
			if ($#argv == 0) goto arg1err;
			set BOLD_stem = $argv[1]; shift;
			breaksw
			
		#skip
		case "-skip":
			if ($#argv == 0) goto arg1err;
			set skip = $argv[1]; shift;
			breaksw

		#echo number
		case "-echo_number":
			if ($#argv == 0) goto arg1err;
			set echo_number = $argv[1]; shift;
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

if ( $skip < 0 ) then
	echo "ERROR: Can not skip negative frame"
	exit 1;
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
	CBIG_preproc_skip.csh
	
DESCRIPTION:
	This function removes the first several frames of fMRI data. For multi-echo data, 
	we will skip the first <num_frames> frames for each echo.

REQUIRED ARGUMENTS:
	-s  <subject_id>        : name of the subject, e.g. Sub0033_Ses1
	-d  <subject_dir>       : full path to <subject_id>, all preprocessed results 
	                          are stored in <subject_dir>/<subject_id>
	-bld  <bold_runs>       : bold numbers of this subject specified by several 
	                          three digits numbers. If there are more than one bold 
	                          number, use a space as delimiter, e.g. '002 003'. 
	                          NOTE: quote sign is necessary.
	-BOLD_stem  <BOLD_stem> : stem of input file, e.g. if the input file name is
	                          Sub0001_Ses1_bld002_rest.nii.gz, the BOLD_stem will be _rest.
	                          This input file should be stored in 
	                          <subject_dir>/<subject_id>/bold/<run_number>/

	                         
OPTIONAL ARGUMENTS:
	-skip  <num_frames>     : skip first several frames. Default is 4.
	-force                  : update results, if exist then overwrite
	-echo_number <echo_number> : number of echoes. For single echo data, default is 1.
	-help                   : help
	-version                : version

Example:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_skip.csh 
	-s Sub0033_Ses1 -d ~/storage/fMRI_preprocess -bld '002 003' -BOLD_stem _rest -skip 4

	For multi-echo case:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_skip.csh 
	-s Sub005 -d ~/storage/fMRI_preprocess -bld '001' -BOLD_stem _rest 
	-skip 4 -echo_number 3
