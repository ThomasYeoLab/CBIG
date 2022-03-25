#!/bin/csh -f

# Example: 
# csh CBIG_preproc_multiecho_denoise.csh -s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -e 12 30.11 48.22 
#
# This function uses tedana (https://tedana.readthedocs.io/en/stable/index.html) to perform multiecho. 
# 'Tedana' is an ICA-based denoising pipeline built for multi-echo data. 
# Rather than analyzing single-echo time series separately, 'tedana' first optimally combines multiple echoes
# then runs denoising with an ICA-based method.
#
# Written by Lyu Xingyu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set subject = ""       # subject ID
set sub_dir = ""       # directory to subjects
set bold = ""          # bold numbers, e.g. '002 003'
set BOLD_stem = ""     # BOLD stem, e.g. _rest_skip5_stc_mc_residc
set echo_number = ""   # number of total echos
set echo_path = ""     # path of every echo in order
set echo_time = ""     # echo time for each echo in order
set nocleanup = 0      # default clean up intermediate files
set lowmem = 0         # low-memory processing. May increase workflow duration. Default is false
set multiecho_stem = "_me"
set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

# Print help and version
set VERSION = '$Id: CBIG_preproc_multiecho_denoise.csh v 1.0 2021/10/29'

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

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

cd $sub_dir/$subject

###########################
# Create log file
###########################
if (! -e logs) then
	mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_multiecho_denoise.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[Multiecho]: logfile = $LF"
echo "multiecho" >> $LF
echo "[CMD]: CBIG_preproc_multiecho_denoise.csh $cmdline"   >>$LF

##########################
# Specify BOLD folder, it contains each run folder (like 002 003)
##########################
set boldfolder = "$sub_dir/$subject/bold"
echo "[ME]: boldfolder = $boldfolder" |& tee -a $LF

cd $boldfolder

##########################
# Optimally combined echoes and perform multi-echo ICA (ME-ICA) (TEDANA)
##########################
echo "=====================combine echos and perform multi-echo ICA using tedana ======================" |& tee -a $LF
foreach runfolder ($bold)
	echo ">>> Run: $runfolder"
	pushd $runfolder
	 
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	set output = "${BOLD}${multiecho_stem}.nii.gz"
	if ( -e $output) then
		echo "[ME]: $output already exists." |& tee -a $LF
	else
		set echo_times=($echo_time:as/,/ /)
		set i = 1
		set echo_paths = ""
		
		while ($i<= $echo_number)
			set echo_paths = ($echo_paths "$subject"_bld$runfolder"_e$i"$BOLD_stem".nii.gz")
			@ i++
		end
		
		if ( $#echo_paths != "$echo_number" || $#echo_times != "$echo_number") then
			echo "ERROR: inconsistency in number of echo times or paths with echo number" |& tee -a $LF
			exit 1;
		endif
		if ( -e ME_intermediate) then
			rm -rf ME_intermediate
		endif
		mkdir -p ME_intermediate
		set CMD="tedana -d $echo_paths -e $echo_times --out-dir $boldfolder/$runfolder/ME_intermediate"
		if ( $lowmem == 1 ) then 
			set CMD = "$CMD --lowmem"
		endif
		echo $CMD |& tee -a $LF
		eval $CMD >> $LF

		# We keep images with and without MEICA for QC purpose. And remove other intermediate files if not specified.
		rsync ME_intermediate/desc-optcomDenoised_bold.nii.gz $output
		rsync ME_intermediate/desc-optcomDenoised_bold.nii.gz \
		$sub_dir/$subject/qc/bld${runfolder}_desc-optcomDenoised_bold.nii.gz
		rsync ME_intermediate/desc-optcom_bold.nii.gz $sub_dir/$subject/qc/bld${runfolder}_desc-optcom_bold.nii.gz
		if ( $nocleanup != 1 ) then 
			rm -rf ME_intermediate
		endif
		
	endif
	popd
end
echo "====================== Multi-echo ICA finished. ======================" |& tee -a $LF

#########################
# Output last commit of current function 
#########################
# Check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	pushd ${CBIG_CODE_DIR}
	git log -1 -- ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/\
	CBIG_preproc_multiecho_denoise.csh >> $LF
	popd
endif

echo "******************************************************************************" |& tee -a $LF
echo ""

exit 0;

###########################
##======Pass the arguments
###########################
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
			
		# bold number, e.g. '002 003'
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set bold = ($argv[1]); shift;
			breaksw

		#BOLD stem
		case "-BOLD_stem":
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = $argv[1]; shift;
			breaksw
			
		# echo number
		case "-echo_number":
			if ( $#argv == 0 ) goto arg1err;
			set echo_number = $argv[1]; shift;
			breaksw
			
		# echo time
		case "-echo_time":
			if ( $#argv == 0 ) goto arg1err;
			set echo_time = $argv[1]; shift;
			breaksw	
		
		#cleanup intermediate files
		case "-nocleanup":
			set nocleanup = 1;
			breaksw

		#cleanup intermediate files
		case "-lowmem":
			set lowmem = 1;
			breaksw

		default:
			echo ERROR: Flag $flag unrecognized.
			echo $cmdline
			exit 1
			breaksw
	endsw
end

goto parse_args_return;

#############################
##======check passed parameters
#############################
check_params:

if ( $#subject == 0 ) then
	echo "ERROR: subject not specified"
	exit 1;
endif
 
if ( $#sub_dir == 0 ) then
	echo "ERROR: path to subject folder not specified"
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

if ( $#echo_time == 0 ) then
	echo "ERROR: echo time not specified"
	exit 1;
endif

if ( $#echo_number == 0 ) then
	echo "ERROR: echo number not specified"
	exit 1;
endif
			
goto check_params_return;

##############################			
##======Error message
##############################
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
	CBIG_preproc_multiecho_denoise.csh

DESCRIPTION:
	This function uses tedana to conduct multiecho preprocessing. 

REQUIRED ARGUMENTS:
	-s            <subject>     : fMRI subject id, enter each echo in order and separate with ","! space will incur error!
	-d            <sub_dir>     : absolute path to <subject>, all preprocessed results 
	                              are stored in <sub_dir>/<subject>
	-echo_number  <echo_number> : number of echoes
	-echo_time    <echo_time>   : comma separated echo time (TE) of each echo (in miliseconds). 
								  must be in the same order as input echo path
	-bld          <bold>        : bold numbers of this subject specified by several three 
	                              digits numbers. If there are more than one bold number, 
	                              use a space as delimiter, e.g. '002 003'. 
	                              NOTE: quote sign is necessary.
	-BOLD_stem    <BOLD_stem>   : stem of input file, e.g. if the input file name is
	                              Sub0001_Ses1_bld002_rest_skip5_stc_mc.nii.gz, the BOLD_stem will be _rest_skip5_stc_mc.
	                              This input file should be stored in 
	                              <sub_dir>/<subject>/bold/<run_number>/
OPTIONAL ARGUMENTS:
	-help                       : help
	-version                    : version
	-lowmem                     : Select this flag to enable TEDANA low-memory processing, 
								  which will save ~15% of the memory usage. May increase workflow duration. 
	-nocleanup                  : Select this flag to keep all intermediate files. By default, we will only keep the 
	                              optimally-combined and optimally-combined+MEICA volumes for QC purpose, 
	                              and remove all other intermediate files. 
	                              For a complete list of tedana output, please refer to this link: 
	                              https://tedana.readthedocs.io/en/stable/outputs.html

OUTPUTS:
	This function will output NIFTI volumes <subject>_bld<run_number><BOLD_stem>_me.nii.gz in folder 
	<sub_dir>/<subject>/bold/<run_number>.

Example:
	csh CBIG_preproc_multiecho_denoise.csh -s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -bld
	'002 003' -BOLD_stem _rest_skip5_stc_mc_sdc -echo_number 3 -echo_time 12,30.11,48.22
