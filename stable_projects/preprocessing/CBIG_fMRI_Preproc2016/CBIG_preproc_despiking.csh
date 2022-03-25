#!/bin/csh -f

# Example: 
# csh CBIG_preproc_despiking.csh -s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -bld '002 003' \
# -BOLD_stem _rest_skip5_stc_mc_residc 
#
# This function uses AFNI 3dDespike to perform despiking. 
#
# Written by Xing Qian and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#BOLD: basename of each run input
#boldfolder: directory of /bold
#bold: all runs under /bold folder

set subject = ""       # subject ID
set sub_dir = ""       # directory to subjects
set bold = ""          # bold numbers, e.g. '002 003'
set BOLD_stem = ""     # BOLD stem, e.g. _rest_skip5_stc_mc_residc
set dspk_suffix = "_dspk"

# Print help and version
set VERSION = '$Id: CBIG_preproc_despiking.csh v 1.0 2017/09/22'

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
# create log file
###########################
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_despiking.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[DSPK]: logfile = $LF"
echo "Despiking" >> $LF
echo "[CMD]: CBIG_preproc_despiking.csh $cmdline"   >>$LF

##########################
# specify BOLD folder, it contains each run folder (like 002 003)
##########################
set boldfolder = "$sub_dir/$subject/bold"
echo "[DSPK]: boldfolder = $boldfolder" |& tee -a $LF

pushd $boldfolder

echo "===================== Despiking, using AFNI ======================" |& tee -a $LF

foreach runfolder ($bold)
	echo ">>> Run: $runfolder"
	pushd $runfolder
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	set output = "${BOLD}${dspk_suffix}.nii.gz"
	if(-e $output) then
		echo "[DSPK]: $output already exists." |& tee -a $LF
	else
#########################
# despiking now!
#########################
		set cmd = (3dDespike -overwrite -prefix ${output} ${BOLD}.nii.gz)
		echo $cmd |& tee -a $LF
		eval $cmd >> $LF
	
	endif
	popd
end

popd
echo "====================== Despiking finished. ======================" |& tee -a $LF

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	pushd ${CBIG_CODE_DIR}
	git log -1 -- ${CBIG_CODE_DIR}stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_despiking.csh >> $LF
	popd
endif

echo "******************************************************************************"
echo ""

exit 0;

###########################
##======pass the arguments
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
			set BOLD_stem = "$argv[1]"; shift;
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
# usage exit
#########################################
BEGINHELP

NAME:
	CBIG_preproc_despiking.csh

DESCRIPTION:
	This function uses AFNI 3dDespike to conduct despiking.

REQUIRED ARGUMENTS:
	-s            <subject>   : fMRI subject id
	-d            <sub_dir>   : absolute path to <subject>, all preprocessed results 
	                            are stored in <sub_dir>/<subject>
	-bld          <bold>      : bold numbers of this subject specified by several three 
	                            digits numbers. If there are more than one bold number, 
	                            use a space as delimiter, e.g. '002 003'. 
	                            NOTE: quote sign is necessary.
	-BOLD_stem    <BOLD_stem> : stem of input file, e.g. if the input file name is
	                            Sub0001_Ses1_bld002_rest_skip5_stc_mc_residc.nii.gz, 
								the BOLD_stem will be _rest_skip5_stc_mc_residc.
	                            This input file should be stored in 
	                            <sub_dir>/<subject>/bold/<run_number>/

OUTPUTS:
	This function will output NIFTI volumes <subject>_bld<run_number><BOLD_stem>_dspk.nii.gz in folder 
	<sub_dir>/<subject>/bold/<run_number>.

Example:
	csh CBIG_preproc_despiking.csh -s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -bld
	'002 003' -BOLD_stem _rest_skip5_stc_mc_residc 



