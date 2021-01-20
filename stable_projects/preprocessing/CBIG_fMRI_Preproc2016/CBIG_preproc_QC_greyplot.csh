#! /bin/csh -f

# Create grey plot (plot grey matter timeseries).
#
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_QC_greyplot.csh, v 1.0 2016/06/09 $'

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

set subject = ""           # subject ID
set sub_dir = ""           # directory to subjects
set anat_s = ""            # recon-all folder
set anat_dir = ""          # directory to recon-all folder
set bold = ""              # bold number, like '002 003'
set BOLD_stem = ""         # BOLD stem, e.g. _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08
set reg_stem = ""          # registration stem, e.g. _rest_skip4_stc_mc_reg
set mc_stem = ""           # the stem of input file of motion correction step, e.g. _rest_skip4_stc
set grey_vox_fac = 200     # a factor used to adjust the height of grey matter timeseries subplot
set tp_fac = 0.3           # a factor used to adjust the width of grey matter timeseries subplot
set FD_th = 0.2            # the threshold of FD in censoring
set DV_th = 50             # the threshold of DV in censoring

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

###############################
# check if matlab exists
###############################
set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
if ($status) then
	echo "ERROR: could not find matlab"
	exit 1;
endif

cd $sub_dir/$subject


###############################
# Create log file
###############################
set mc = $sub_dir/$subject"/bold/mc"
set qc = $sub_dir/$subject"/qc"
if (! -e $mc) then
     mkdir -p $mc
endif
if (! -e $qc) then
     mkdir -p $qc
endif
if (! -e logs) then
     mkdir -p logs
endif
echo 6
set LF = $sub_dir/$subject/logs/CBIG_preproc_QC_greyplot.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[Greyplot]: logfile = $LF"
echo "Greyplot" >> $LF
echo "[CMD]: CBIG_preproc_QC_greyplot.csh $cmdline"   >>$LF

set boldfolder = "$sub_dir/$subject/bold"
echo "[Greyplot]: boldfolder = $boldfolder" |& tee -a $LF
echo "[Greyplot]: bold = $bold" |& tee -a $LF


###############################
# Create whole brain & grey matter mask, if not exist
###############################
echo "========================= Create whole brain & Grey Matter Mask =========================" |& tee -a $LF
if( -e $boldfolder/mask/${subject}.func.gm.nii.gz && -e $boldfolder/mask/${subject}.brainmask.bin.nii.gz ) then
	# if grey matter mask and whole brain mask already exist, do nothing
	echo "[Greyplot]: Grey mask and whole brain mask of subject $subject already exist." |& tee -a $LF
else
	# if one or both of GM and whole brain masks does not exist, create them.
	set cmd = "${root_dir}/CBIG_preproc_create_mask.csh -s $subject -d $sub_dir -anat_s ${anat_s} -anat_d ${anat_dir}"
	set cmd = "$cmd -bld '${bold}' -REG_stem ${reg_stem} -MASK_stem ${BOLD_stem} -whole_brain -gm"
	echo $cmd |& tee -a $LF
	eval $cmd
endif
echo "===================== Creating whole brain & Grey Matter Mask finished ===================" |& tee -a $LF
echo "" |& tee -a $LF


###############################
# generate global signal text file
###############################
echo "================================ Generate global signal text file ==============================" |& tee -a $LF
set regress_folder = "$boldfolder/regression"
set ROI_regressors_list = "$regress_folder/ROI_regressors_list.txt"
# if nuisance regression step was performed, global signal should be computed on the volume before nuisance regression
if( -e $ROI_regressors_list ) then   
	set regression_done = 1    
	echo "[Greyplot]: Nuisance regression was one of previous preprocessing steps." |& tee -a $LF
	echo "[Greyplot]: Global signal in the plot will be generated from the input volume of regression step." |& tee -a $LF
	
	foreach curr_run ($bold)
		if ( -e $qc/$subject"_bld"${curr_run}_WB.txt ) then
			rm $qc/$subject"_bld"${curr_run}_WB.txt
		endif
	end
	
	set ROI_regressors_files = (`cat $ROI_regressors_list`)
	set bold_list = ($bold)
	set first_run = `echo $bold_list[1]`
	set first_file = `echo $ROI_regressors_files[1]`
	# check if global signal was included in the nuisance regression
	set prefix = "$regress_folder/$subject"_bld"${first_run}_WB"
	set status1 = 0;
	set status2 = 0;
	set status3 = 0;
	set status4 = 0;
	if( "$first_file" == "${prefix}_regressors.txt" ) then 
		set status1 = 1; 
	endif
	if( "$first_file" == "${prefix}_WM_regressors.txt" ) then 
		set status2 = 1; 
	endif
	if( "$first_file" == "${prefix}_CSF_regressors.txt" ) then 
		set status3 = 1; 
	endif
	if( "$first_file" == "${prefix}_WM_CSF_regressors.txt" ) then
		set status4 = 1; 
	endif
	if( $status1 || $status2 || $status3 || $status4 ) then   # GS was one of the regressor
		echo "[Greyplot]: global signal was generated in regression step. Extract it for later plotting." |& tee -a $LF
		foreach curr_ROI_file ($ROI_regressors_files)
			echo $curr_ROI_file
			set curr_run = `echo $curr_ROI_file | xargs basename | sed -e "s/^${subject}_bld//" | awk '{print substr($0, 0, 3)}'`
			echo $qc/$subject"_bld"${curr_run}_WB.txt
			awk -F " " '{print $1}' $curr_ROI_file >> $qc/$subject"_bld"${curr_run}_WB.txt
		end
	else                                                       # GS was not in the regressor, generate now
		echo "[Greyplot]: global signal (GS) was NOT included in regression step. Generate GS text now." |& tee -a $LF
		set wb_mask = "$boldfolder/mask/$subject.brainmask.bin.nii.gz"
		set GS_list = "$qc/GS_list.txt"
		if( -e $GS_list ) then
			rm $GS_list
		endif
		touch $GS_list
		foreach curr_run ($bold)
			echo $qc/$subject"_bld"${curr_run}_WB.txt >> $GS_list
		end
		
		set fMRI_list = "$regress_folder/fMRI_list.txt"
		set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r )
		set cmd = ($cmd '"' 'addpath(fullfile('"'"$root_dir"'"','"'"utilities"'"'))'; )
		set cmd = ($cmd CBIG_preproc_create_ROI_regressors "'"$fMRI_list"'" "'"$GS_list"'" )
		set cmd = ($cmd "'"$wb_mask"'" "'""'" "'""'" "'"0"'"; exit '"' );
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		rm $GS_list
	endif
# if nuisance regression step was not performed, global signal can be computed on the final volume
else 
	set regression_done = 0
	echo "[Greyplot]: Nuisance regression was not one of previous preprocessing steps." |& tee -a $LF
	echo "[Greyplot]: Global signal in the plot will be generated based on the volume from last preprocessing step." |& tee -a $LF
endif


###############################
# call matlab function to plot
###############################
echo "================================ Create greyplot ==============================" |& tee -a $LF
set wb_mask = "$boldfolder/mask/$subject.brainmask.bin.nii.gz"
set gm_mask = $boldfolder/mask/$subject.func.gm.nii.gz
foreach runfolder ($bold)
	set fmri_file = $boldfolder/$runfolder/${subject}_bld${runfolder}${BOLD_stem}.nii.gz
	set FD_file = $mc/${subject}_bld${runfolder}${mc_stem}_motion_outliers_FDRMS
	set DV_file = $mc/${subject}_bld${runfolder}${mc_stem}_motion_outliers_DVARS
	set output = $qc/${subject}_bld${runfolder}${BOLD_stem}_greyplot.png
	
	if ( "$regression_done" == 0 ) then
		set cmd = ( $MATLAB -nodesktop  -nosplash -r '"' 'addpath(genpath('"'"${root_dir}'/utilities'"'"'))'; )
		set cmd = ($cmd CBIG_preproc_QC_greyplot "'"$fmri_file"'"  "'"$FD_file"'"  "'"$DV_file"'"  "'"$output"'" )  
		set cmd = ($cmd "'"GM_mask"'" "'"$gm_mask"'"  "'"WB_mask"'" "'"$wb_mask"'" "'"grey_vox_factor"'" "'"$grey_vox_fac"'")
		set cmd = ($cmd "'"tp_factor"'" "'"$tp_fac"'" "'"FD_thres"'" "'"$FD_th"'" "'"DV_thres"'" "'"$DV_th"'"; exit; '"')
	else
		set GS_txt = $qc/$subject"_bld"${runfolder}_WB.txt
		
		set cmd = ( $MATLAB -nodesktop  -nosplash -r '"' 'addpath(genpath('"'"${root_dir}'/utilities'"'"'))'; )
		set cmd = ($cmd CBIG_preproc_QC_greyplot "'"$fmri_file"'"  "'"$FD_file"'"  "'"$DV_file"'"  "'"$output"'" )  
		set cmd = ($cmd "'"GM_mask"'" "'"$gm_mask"'"  "'"WBS_text"'" "'"$GS_txt"'" "'"grey_vox_factor"'" "'"$grey_vox_fac"'")
		set cmd = ($cmd "'"tp_factor"'" "'"$tp_fac"'" "'"FD_thres"'" "'"$FD_th"'" "'"DV_thres"'" "'"$DV_th"'"; exit; '"')
	endif
	echo $cmd |& tee -a $LF
	eval $cmd |& tee -a $LF
end
echo "=========================== Creatting greyplot finished =========================" |& tee -a $LF
echo "" |& tee -a $LF



#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	pushd ${CBIG_CODE_DIR}
	git log -1 -- ${CBIG_CODE_DIR}/table_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_QC_greyplot.csh >> $LF
	popd
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
			
		case "-MC_stem":
			if ( $#argv == 0 ) goto arg1err;
			set mc_stem = "$argv[1]"; shift;
			breaksw
			
		case "-grey_vox_fac"
			if ( $#argv == 0 ) goto arg1err;
			set grey_vox_fac = "$argv[1]"; shift;
			breaksw
			
		case "-tp_fac"
			if ( $#argv == 0 ) goto arg1err;
			set tp_fac = "$argv[1]"; shift;
			breaksw
			
		case "-FD_th"
			if ( $#argv == 0 ) goto arg1err;
			set FD_th = "$argv[1]"; shift;
			breaksw
			
		case "-DV_th"
			if ( $#argv == 0 ) goto arg1err;
			set DV_th = "$argv[1]"; shift;
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
	echo "ERROR: reg stem not specified"
	exit 1;
endif

if ( "$mc_stem" == "" ) then
	echo "ERROR: mc stem not specified"
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
	CBIG_preproc_QC_greyplot.csh

DESCRIPTION:
	This function
	  1) Create whole brain and grey matter mask, if not exist.
	  2) Call matlab function to create the greyplot.
	
	Greyplot contains 4 subplots: (1) framewise displacement trace;
	(2) DVARS trace; (3) global signal; and (4) grey matter voxels 
	timseries.
	
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
	                              Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08.nii.gz,  
	                              then <BOLD_stem> = _rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08. 
	                              This input file is assumed to be stored in 
	                              <sub_dir>/<subject>/bold/<run_number>.
	-REG_stem      reg_stem     : stem of T1-T2* registration file. E.g. if the registration 
	                              file is Sub0001_Ses1_bld002_rest_skip4_stc_mc_reg.dat, then 
	                              <reg_stem> = _rest_skip4_stc_mc_reg. The registration file is 
	                              assumed to be stored in <sub_dir>/<subject>/bold/<run_number>.
	-MC_stem       mc_stem      : stem of input volume of motion correction step. E.g. if the 
	                              input file of motion correction step is 
	                              Sub0001_Ses1_bld002_rest_skip4_stc.nii.gz, then <mc_stem> = 
	                              _rest_skip4_stc. This file is assumed to be stored in 
	                              <sub_dir>/<subject>/bold/<run_number>.

OPTIONAL ARGUMENTS:
	-grey_vox_fac  grey_cox_fac : it is a factor used to adjust the height of grey matter 
	                              timeseries subplot. The height is defined as: 
	                              (number of grey matter voxels) / (grey_vox_fac).
	                              Default grey_vox_fac is 200.
	                              
	-tp_fac        tp_fac       : it is a factor used to adjust the width of grey matter 
	                              timeseries subplot. The width is defined as:
	                              (number of timepoints) / (tp_fac).
	                              Default tp_fac is 0.3.
	-FD_th         FD_th        : it is the threshold for FD used in detecting censored frames. 
	                              Default is 0.2.
	-DV_th         DV_th        : it is the threshold for DVARS used in detecting censored frames.
	                              Default is 50.
	-help                       : help
	-version                    : version

OUTPUTS:
	For each run, the greyplot is saved as
	<subject_dir>/<subject_id>/qc/<subject>_bld<run_number><BOLD_stem>_greyplot.png

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_QC_greyplot.csh 
	-s Sub0001_Ses1 -d /storage/FMRI_preprocess -bld '002 003' -BOLD_stem 
	_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08 -REG_stem _rest_skip4_stc_mc_reg
	-MC_stem _rest_skip4_stc -grey_vox_fac 200 -tp_fac 0.3


