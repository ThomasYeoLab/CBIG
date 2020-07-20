#!/bin/csh -f

# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslslicetimer.csh 
#	-s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4 -slice_order
#	${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/example_slice_order.txt
#
# This function uses FSL slicetimer to conduct slice timing correction. 
# One slice timing file or slice order file is needed. 
# If no one is passed in, the programme will generate a slice order file automatically. The default one is Siemens.
#
# Written by Jingwei Li.
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#BOLD: basename of each run input
#boldfolder: directory of /bold
#bold: all runs under /bold folder

set subject = ""       # subject ID
set sub_dir = ""       # directory to subjects
set bold = ""          # bold numbers, e.g. '002 003'
set BOLD_stem = ""     # BOLD stem, e.g. _rest
set stc_suffix = "_stc"

# Print help and version
set VERSION = '$Id: CBIG_preproc_fslslicetimer.csh v 1.0 2016/05/26'

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

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:


set currdir = `pwd`
cd $sub_dir/$subject



###########################
# create log file
###########################
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_fslslicetimer.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[STC]: logfile = $LF"
echo "Slice-Time Correction" >> $LF
echo "[CMD]: CBIG_preproc_fslslicetimer.csh $cmdline"   >>$LF

##########################
# specify BOLD folder, it contains each run folder (like 002 003)
##########################
set boldfolder = "$sub_dir/$subject/bold"
echo "[STC]: boldfolder = $boldfolder" |& tee -a $LF

pushd $boldfolder

echo "===================== Slice time correction, using fsl ======================" |& tee -a $LF

set run_count = 0
foreach runfolder ($bold)
	echo ">>> Run: $runfolder"
	set run_count = `expr $run_count + 1`
	pushd $runfolder
	set BOLD = ${subject}"_bld${runfolder}${BOLD_stem}"
	set output = "${BOLD}_stc.nii.gz"
	if(-e $output) then
		echo "[STC]: $output already exists." |& tee -a $LF
	else
		set nslices = `fslval $BOLD.nii.gz dim3`
		
		#########################
		# if slice acquisition direction is not passed in, create a default one (z=3 from superior to inferior or vice versa)
		#########################
		if(! $?direction ) then
			set direction = 3; 
		endif
		
		##########################
		# if slice-timing file in passed in, check how many columns (each column inidicates slice-timing of one run) are there
		# if > 1, split multiple columns into multiple files
		# slice timings are normalized to TR
		##########################
		set auto_st_file_flag = 0;
		if( $?st_file_in ) then
			set MRI_info = `mri_info ${BOLD}.nii.gz`
			# extract TR from header information
			set TR = `echo $MRI_info | grep -o 'TR: \(.*\)' | awk -F " " '{print $2}'`
			if ( "$TR" == "" ) then
				echo "ERROR: TR cannot be extracted from image header information" |& tee -a $LF
				exit 1
			endif
			# check the unit for TR is msec or sec
			set unit = `echo $MRI_info | grep -o 'TR: \(.*\)' | awk -F " " '{print $3}'`
			set unit = `echo $unit | awk -F "," '{print $1}'`
			# convert the unit of TR to sec (required by FSL slicetimer)
			if ( "$unit" == "msec" ) then
				set TR = `echo "scale=4; $TR/1000" | bc`
			else
				echo "ERROR: Unknown TR unit" |& tee -a $LF
				exit 1
			endif
			set num_col = `awk -F ' ' '{print NF; exit}' $st_file_in`
			set num_run = `echo $#bold`
			if( $num_col == 1 ) then
				set st_file = $st_file_in
				set st_fraction_file = $boldfolder/$runfolder/tmp_stc/tmp_st_fraction.txt
				mkdir -p tmp_stc
				foreach slice_timing ("`cat $st_file`")
					echo "scale=4; 0.5 - $slice_timing/$TR" | bc >> $st_fraction_file
				end	
				set st_file = $st_fraction_file	
			else if ( $num_col == $num_run ) then
				set st_file = $boldfolder/$runfolder/tmp_stc/tmp_st.txt
				mkdir -p tmp_stc
				cut -f $run_count -d ' ' $st_file_in > $st_file
				set st_fraction_file = $boldfolder/$runfolder/tmp_stc/tmp_st_fraction.txt
				foreach slice_timing ("`cat $st_file`")
					echo "scale=4; 0.5 - $slice_timing/$TR" | bc >> $st_fraction_file
				end	
				set st_file = $st_fraction_file	
				set auto_st_file_flag = 1
			else
				echo "ERROR: Number of columns in slice-timing file is not consistent with total number of runs." \
				|& tee -a $LF
				exit 1
			endif
		endif

	
		#########################
		# if neither slice order file nor slice timing file is passed in, 
		# create a default slice order file (Siemens)
		#########################
		if( ! $?st_file && ! $?so_file ) then
			set auto_so_file_flag = 1;
			echo "  WARNNING: Slice order file not specified, create a temporary one" |& tee -a $LF
			mkdir -p tmp_stc
			set so_file = $boldfolder/$runfolder/tmp_stc/tmp_so.txt
			if( -e ${so_file} ) then
				rm $so_file
			endif
			
			set n = `expr $nslices % 2`
			if ( $n == 1 ) then
				@ nthslice = 1
				while($nthslice <= $nslices)
					echo $nthslice >> $so_file
					@ nthslice = $nthslice + 2;
				end
				@ nthslice = 2
				while($nthslice <= $nslices)
					echo $nthslice >> $so_file
					@ nthslice = $nthslice + 2;
				end
			else
				@ nthslice = 2
				while($nthslice <= $nslices)
					echo $nthslice >> $so_file
					@ nthslice = $nthslice + 2;
				end
				@ nthslice = 1
				while($nthslice <= $nslices)
					echo $nthslice >> $so_file
					@ nthslice = $nthslice + 2;
				end
			endif
		else
			set auto_so_file_flag = 0;
		endif
		
		
		#########################
		# slice timing correction now!
		#########################
		set cmd = (slicetimer -i ${BOLD}.nii.gz -o ${BOLD}${stc_suffix}.nii.gz -d $direction)
		
		if( $?st_file ) then
			echo "  ------------------ Slice Timing (in unit of TR) -------------------" |& tee -a $LF
			cat $st_file | tr '\n' ' ' |& tee -a $LF
			echo "" |& tee -a $LF
			echo "  -------------------------------------------------------------------" |& tee -a $LF
			
			set cmd = ($cmd --tcustom=$st_file)
			echo $cmd |& tee -a $LF
			eval $cmd >> $LF
		else
			echo "  ---------------------- Slice Order --------------------------------" |& tee -a $LF
			cat $so_file | tr '\n' ' ' |& tee -a $LF
			echo "" |& tee -a $LF
			echo "  -------------------------------------------------------------------" |& tee -a $LF
			
			set cmd = ($cmd --ocustom=$so_file)
			echo $cmd |& tee -a $LF
			eval $cmd >> $LF
		endif
		
		#########################
		# remove the automatically generated default slice order file, if there is one
		# remove the split slice timing file of current run 
		# if input slice-timing file contains multiple columns of slice timing
		#########################
		if ( $auto_so_file_flag == 1 || $auto_st_file_flag == 1 ) then
			rm -r tmp_stc
			unset so_file
			unset st_file
		endif
	endif
	popd
end

popd
echo "====================== Slice time correction finished. ======================" |& tee -a $LF

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslslicetimer.csh\
>> $LF
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
		
		# slice aquisition direction
		case "-direction":
			if ( $#argv == 0 ) goto arg1err;
			set direction = $argv[1]; shift;
			breaksw
	
		# slice order file (each row is one number)
		case "-slice_order":
			if ( $#argv == 0 ) goto arg1err;
			set so_file = $argv[1]; shift;
			breaksw
			
		# slice timing file
		case "-slice_timing":
			if ( $#argv == 0 ) goto arg1err;
			set st_file_in = $argv[1]; shift;
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
	CBIG_preproc_fslslicetimer.csh

DESCRIPTION:
	This function uses FSL slicetimer to conduct slice timing correction. One slice order file is needed. 
	If the slice order file is not passed in, the programme will generate one automatically. The default 
	one is the interleaved order (Siemens: if the number of slices is odd, the ordering is 1, 3, 5, ..., 
	2, 4, 6, ...; if the number of slices is even, the ordering is 2, 4, 6, ..., 1, 3, 5, ...). The output 
	stem is BOLD stem plus '_stc'

REQUIRED ARGUMENTS:
	-s            <subject>   : fMRI subject id
	-d            <sub_dir>   : absolute path to <subject>, all preprocessed results 
	                            are stored in <sub_dir>/<subject>
	-bld          <bold>      : bold numbers of this subject specified by several three 
	                            digits numbers. If there are more than one bold number, 
	                            use a space as delimiter, e.g. '002 003'. 
	                            NOTE: quote sign is necessary.
	-BOLD_stem    <BOLD_stem> : stem of input file, e.g. if the input file name is
	                            Sub0001_Ses1_bld002_rest_skip4.nii.gz, the BOLD_stem will be _rest_skip4.
	                            This input file should be stored in 
	                            <sub_dir>/<subject>/bold/<run_number>/

OPTIONAL ARGUMENTS:
	-direction    <direction> : slice aquisition direction. 1 means x axis representing Right-Left, 
	                            2 means y axis representing Anterior-Posterior, 
	                            3 means z axis representing Superior-Inferior.
	-slice_timing <st_file>   : slice timing file (absolute path). It can contain one column or multiple columns 
	                            of slice timing (separated by a space, see example_slice_timing.txt within the same 
	                            folder). 
	                            If there is only one column of slice timing, it is assumed that every run follows 
	                            the same slice timing;  
	                            if there are multiple columns of slice timing, each column corresponds to 
	                            one of the runs given by -bld option.
	                            N-th row corresponds to the time point when N-th slice was sampled.
	                            NOTE: If both -slice_timing and -slice_order options are not used, 
	                            this script will generate a default slice order file following Siemens slice order 
	                            (see description of -slice_order option).
	-slice_order  <so_file>   : slice order file (absolute path), where each row is one number. 
	                            If the user did not pass in a slice order file or a slice timing file,  
	                            this function will create a default one where the slice order is 
	                            following Siemens: if the number of slices is odd, the ordering is 
	                            1, 3, 5, ..., 2, 4, 6, ...; if the number of slices is even, the 
	                            ordering is 2, 4, 6, ..., 1, 3, 5, ....
	                            For the example of slice order file, please see 
	                            ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/\
	                            example_slice_order.txt
OUTPUTS:
	This function will output NIFTI volumes <subject>_bld<run_number><BOLD_stem>_stc.nii.gz in folder 
	<sub_dir>/<subject>/bold/<run_number>.

Example:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslslicetimer.csh 
	-s Sub0001_Ses1 -d ~/storage/fMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4 -slice_order
	${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/example_slice_order.txt


Written by Jingwei Li.
Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

