#! /bin/csh -f

#############################################
# Motion Correction and Outlier Detection
#############################################
#############################################
# In this script, we: 
# 1) use mcflirt to do the motion correction, 
# 2) use fsl_motion_outliers to obtain FDRMS and DVARS
# 3) use DV_th and FD_th as threhold to find outliers, the default DV_th is 50, default of FD_th is 0.2
# 4) frames either above FD_th or above DV_th will be removed
# 5) 1 frame before 2 frames after each removed frame will also be removed, as well as kept segments of data lasting 
#    fewer than $discard_seg contiguous frames will be removed. If not specified, the default is 5.
# 6) discard the run which has more than <rm_run_th>% of frames being removed
# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslmcflirt_outliers.csh 
#	-s Sub0033_Ses1 -d ~/storage/FMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4_stc -nframe 0 
#	-FD_th 0.2 -DVARS 50 -discard-run 50 -rm-seg 5
#############################################
# Author: RU(BY) KONG
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_fslmcflirt_outliers.csh, v 1.0 2016/06/09 $'

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
set sub_dir = ""
set zpdbold = ""
set BOLD_stem = ""
set nframe = 0;  # Perform mcflirt with the n-th frame
set fd_th = 0.2; # Default FDRMS threshold is 0.2
set dv_th = 50;  # Default DVARS threshold is 50
set force = 0;   # Default if file exist, then skip this step.
set nocleanup = 0; # Default clean up intermediate file
set discard_seg = 5; # Default will remove kept segments of data lasting fewer than 5 contiguous frames
set rm_run_th = 50; # Default will discard the run which has more than 50% of frames being removed
set spline_final = 0; # spline_final flag in mcflirt, 0 for trilinear interpolation, 1 for spline interpolation

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

#############################################
# BOLD Information
#############################################

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
set LF = $sub_dir/$subject/logs/CBIG_preproc_fslmcflirt_outliers.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[MC]: logfile = $LF"
echo "Motion Correction" >> $LF
echo "[CMD]: CBIG_preproc_fslmcflirt_outliers.csh $cmdline"   >>$LF

set boldfolder = "$sub_dir/$subject/bold"
echo "[MC]: boldfolder = $boldfolder" |& tee -a $LF
echo "[MC]: zpdbold = $zpdbold" |& tee -a $LF

##base bold
cd $boldfolder
set base_bold = $boldfolder/$zpdbold[1]
echo "[MC]: base_bold = $base_bold" |& tee -a $LF



#############################################
# Generate Template use n-th frame of the first run
#############################################

echo "=======================Generate template.nii.gz..(n-th frame of the base_bold, default is 1st frame)=======================" |& tee -a $LF
pushd $base_bold

set base_boldfile = $subject"_bld"$zpdbold[1]$BOLD_stem
echo "[MC]: use $nframe +1 as the template" |& tee -a $LF
if ( (! -e  $boldfolder/mc_template.nii.gz) || ($force == 1)) then
	#Motion correction template is the nframe = $nframe, default is the first frame
	fslroi $base_boldfile $boldfolder/mc_template $nframe 1 |& tee -a $LF 
else
	echo "[MC]: Template already exists!" |& tee -a $LF
endif	

popd
set template = $boldfolder/mc_template.nii.gz
echo "[MC]: template = $template" |& tee -a $LF
echo "=======================Template generation done!=======================" |& tee -a $LF
echo "" |& tee -a $LF

#############################################
# Merge each run with template
#############################################

echo "=======================Merge each run with template=======================" |& tee -a $LF
foreach curr_bold ($zpdbold)
	pushd $curr_bold
	
	set boldfile = $subject"_bld"$curr_bold$BOLD_stem
	echo "[MC]: boldfile = $boldfile" |& tee -a $LF
	if ( (! -e  $boldfile"_merge.nii.gz") || ($force == 1) ) then
		fslmerge -t $boldfile"_merge" $template $boldfile |& tee -a $LF
	else
		echo "=======================Merge already done!=======================" |& tee -a $LF
	endif
	
	popd
end
echo "=======================Merge done!=======================" |& tee -a $LF
echo "" |& tee -a $LF

#############################################
# Motion correction: mcflirt, split the template from results
#############################################

echo "=======================Motion correction: mcflirt and split=======================" |& tee -a $LF
foreach curr_bold ($zpdbold)
	pushd $curr_bold
	
	set boldfile = $subject"_bld"$curr_bold$BOLD_stem
	if ( (! -e  $boldfile"_mc.nii.gz") || ($force == 1) ) then
		set cmd = "mcflirt -in ${boldfile}_merge.nii.gz -out ${boldfile}_mc -plots -refvol 0 -rmsrel -rmsabs -mats"
		if ($spline_final == 1) then
			set cmd = "$cmd -spline_final"
		endif
		echo $cmd |& tee -a $LF
		eval $cmd >> $LF
		mv $boldfile"_mc.nii.gz" $boldfile"_mc_tmp.nii.gz"
		set numof_tps = `fslnvols $boldfile` 
		fslroi $boldfile"_mc_tmp" $boldfile"_mc" 1 $numof_tps |& tee -a $LF
		rm $boldfile"_mc_tmp.nii.gz"
	else
		echo "=======================mcflirt and split already done!=======================" |& tee -a $LF
	endif
	rm ${boldfile}_mc.mat/MAT_0000
	cat ${boldfile}_mc.mat/MAT* > ${boldfile}_mc.cat
	ln -s $boldfolder/$curr_bold/${boldfile}_mc_rel.rms $qc/${subject}_bld${curr_bold}_mc_rel.rms
	ln -s $boldfolder/$curr_bold/${boldfile}_mc_abs.rms $qc/${subject}_bld${curr_bold}_mc_abs.rms
	ln -s $boldfolder/$curr_bold/${boldfile}_mc_rel_mean.rms $qc/${subject}_bld${curr_bold}_mc_rel_mean.rms
	ln -s $boldfolder/$curr_bold/${boldfile}_mc_abs_mean.rms $qc/${subject}_bld${curr_bold}_mc_abs_mean.rms
	
	popd
end
echo "=======================mcflirt and split done!=======================" |& tee -a $LF
echo "" |& tee -a $LF


#############################################
# Plot mcflirt parameters
#############################################
echo "=========================== Plot mcflirt parameters =============================" |& tee -a $LF
foreach curr_bold ($zpdbold)
	pushd $curr_bold
	
	pwd |& tee -a $LF
	set mc_par_file = "${subject}_bld${curr_bold}${BOLD_stem}_mc.par"
	set mc_abs_rms_file = "${subject}_bld${curr_bold}${BOLD_stem}_mc_abs.rms"
	set mc_rel_rms_file = "${subject}_bld${curr_bold}${BOLD_stem}_mc_rel.rms"
	set outname_prefix = "${subject}_bld${curr_bold}${BOLD_stem}_mc"
	
	set cmd = ( $MATLAB -nodesktop -nodisplay -nosplash -r '"' 'addpath(genpath('"'"${root_dir}'/utilities'"'"'))'; CBIG_preproc_plot_mcflirt_par $mc_par_file $mc_abs_rms_file $mc_rel_rms_file $qc $outname_prefix; exit; '"' );
	echo $cmd |& tee -a $LF
	eval $cmd |& tee -a $LF
	
	popd
end
echo "=========================== Plot mcflirt parameters done =============================" |& tee -a $LF


#############################################
#compute FDRMS and DVARS
#############################################

echo "=======================FSL motion outliers=======================" |& tee -a $LF
foreach curr_bold ($zpdbold)
	pushd $curr_bold
	
	set boldfile = $subject"_bld"$curr_bold$BOLD_stem
	mkdir -p $mc/tmp_outliers/$curr_bold
	
	# Use DVARS as metric
	if ( (! -e $mc/${boldfile}_motion_outliers_DVARS) || ($force == 1) ) then
		echo "[MC]: bold = $curr_bold Perform FSL motion outliers with metric = dvars" |& tee -a $LF
		set cmd = "fsl_motion_outliers -i ${boldfile}_mc -o $mc/${boldfile}_motion_outliers_confound_DVARS -s $mc/${boldfile}_motion_outliers_DVARS -p $mc/${boldfile}_motion_outliers_DVARS -t $mc/tmp_outliers/$curr_bold --dvars --nomoco"
		echo $cmd |& tee -a $LF
		eval $cmd >> $LF
	else
		echo "[MC]: bold = $curr_bold Perform FSL motion outliers with metric = dvars already done!"
	endif
	
	#Use FDRMS as metric
	if ( (! -e $mc/${boldfile}_motion_outliers_FDRMS) || ($force == 1) ) then
		echo "[MC]: bold = $curr_bold Perform FSL motion outliers with metric = fdrms" |& tee -a $LF
	
		set cmd = "fsl_motion_outliers -i $boldfile -o $mc/${boldfile}_motion_outliers_confound_FDRMS -s $mc/${boldfile}_motion_outliers_FDRMS -p $mc/${boldfile}_motion_outliers_FDRMS -t $mc/tmp_outliers/$curr_bold --fdrms"
		echo $cmd |& tee -a $LF
		eval $cmd >> $LF
	else
		echo "[MC]: bold = $curr_bold Perform FSL motion outliers with metric = refrms already done!"
	endif
	
	popd
end
echo "=======================FSL motion outliers done!=======================" |& tee -a $LF
echo "" |& tee -a $LF

#############################################
# Outlier Detection
#############################################

echo "=======================Compute correlation between DVARS and FDRMS and detect the motion outliers=======================" |& tee -a $LF
foreach curr_bold ($zpdbold)
	pushd $curr_bold
	
	set boldfile = $subject"_bld"$curr_bold$BOLD_stem
	if ( (! -e "$qc/${subject}_bld${curr_bold}_FDRMS${fd_th}_DVARS${dv_th}_motion_outliers.txt") || ($force == 1) ) then
		set output = "$qc/${subject}_bld${curr_bold}"
		set dvars_file = "$mc/${boldfile}_motion_outliers_DVARS"
		set fd_file = "$mc/${boldfile}_motion_outliers_FDRMS"
		set cmd = ( $MATLAB -nodesktop -nodisplay -nosplash -r '"' 'addpath(genpath('"'"${root_dir}'/utilities'"'"'))'; CBIG_preproc_DVARS_FDRMS_Correlation $dvars_file $fd_file $output; CBIG_preproc_motion_outliers $dvars_file $fd_file $fd_th $dv_th $discard_seg $output; exit; '"' );
		eval $cmd |& tee -a $LF
	else
		echo "[MC]: Motion outliers detection already created!" |& tee -a $LF
	endif
	echo "[MC]: Motion outliers is in $qc/${subject}_bld${curr_bold}_FDRMS${fd_th}_DVARS${dv_th}_motion_outliers.txt" |& tee -a $LF
	
	popd
end

echo "*********************************************************************" |& tee -a $LF

#############################################
# Clean up intermediate files
#############################################
if ( $nocleanup != 1 ) then
	foreach curr_bold ($zpdbold)
		pushd $curr_bold
		
		set boldfile = $subject"_bld"$curr_bold$BOLD_stem
		set mergefile = $boldfile"_merge.nii.gz"
		rm $mergefile
		
		popd
	end
endif


###############################################
# check if the number of outliers exceeds 50% of total number of frames
###############################################
if ($rm_run == 1) then 
	echo "========= check if the number of outliers exceeds $rm_run_th% of total number of frames for each run =========" |& tee -a $LF
	set bold_file = "$sub_dir/$subject/logs/${subject}.bold"
	echo "[MC]: bold_file = $bold_file" |& tee -a $LF
	foreach runfolder ($zpdbold)
		set outlier_file = "$sub_dir/$subject/qc/${subject}_bld${runfolder}_FDRMS${fd_th}_DVARS${dv_th}_motion_outliers.txt"
		echo "[MC]: outlier_file = $outlier_file" |& tee -a $LF
		set num_zeros = `fgrep -o 0 ${outlier_file} | wc -l`
		echo "[MC]: num_zeros = $num_zeros" |& tee -a $LF
		set num_ones = `fgrep -o 1 ${outlier_file} | wc -l`
		echo "[MC]: num_ones = $num_ones" |& tee -a $LF
		@ num_frames = $num_zeros + $num_ones
		echo "[MC]: num_frames = $num_frames" |& tee -a $LF
		set prop = `awk "BEGIN {print $num_zeros * 100 / $num_frames}"`
		echo "[MC]: prop = $prop" |& tee -a $LF
	
		if( `echo "$prop > $rm_run_th" | bc` ) then
			echo "Run $runfolder has more than $rm_run_th% outliers, remove this run from bold file." |& tee -a $LF
			sed -i -e "s/${runfolder}//g" ${bold_file}
			set bold = `cat ${bold_file}`
		else
			echo "Run $runfolder has less than $rm_run_th% outliers, nothing change to this run." |& tee -a $LF
		endif
	end
	echo "====================== $rm_run_th% outliers threshold check finished ======================" |& tee -a $LF
	echo "" |& tee -a $LF
endif

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslmcflirt_outliers.csh >> $LF
endif

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
			
		#path to subject's preprocess output folder
		case "-d":
			if ( $#argv == 0 ) goto arg1err;
			set sub_dir = $argv[1]; shift;
			breaksw
			
		#bold run number	
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set zpdbold = ($argv[1]); shift;
			breaksw
		#input file stem
		case "-BOLD_stem":
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = "$argv[1]"; shift;
			breaksw
			
		#Perform mcflirt with the n-th frame, the default is the 1st frame, nframe = 0			
		case "-nframe":
			if ( $#argv == 0 ) goto arg1err;
			set nframe = "$argv[1]"; shift;
			breaksw
			
		#FDRMS threshold for motion outlier detection	
		case "-FD_th":
			if ( $#argv == 0 ) goto arg1err;
			set fd_th = "$argv[1]"; shift;
			breaksw
			
		#DVARS threshold for motion outlier detection			
		case "-DV_th":
			if ( $#argv == 0 ) goto arg1err;
			set dv_th = "$argv[1]"; shift;
			breaksw
		
		#update results, if exist then do not generate
		case "-force":
			set force = 1;
			breaksw
			
		#clean up intermediate files
		case "-nocleanup":
			set nocleanup =1;
			breaksw
		
		#throw run which has more than 50% frames detected as outliers
		case "-discard-run":
			set rm_run = 1;
			if ( $#argv == 0 ) goto arg1err;
			set rm_run_th = "$argv[1]"; shift;
			breaksw
			
		#remove kept segments of data lasting fewer than $discard_seg contiguous frames
		case "-rm-seg":
			if ( $#argv == 0 ) goto arg1err;
			set discard_seg = "$argv[1]"; shift;
			breaksw
		
		#spline_final flag to use spline interpolation in mcflirt
		case "-spline_final":
			set spline_final = 1;
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
if ( "$sub_dir" == "" ) then
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
goto check_params_return;

##########################################
# ERROR message
##########################################
	
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

	
exit 0;

#-------- Everything below is printed as part of help --------#
BEGINHELP

NAME: 
	CBIG_preproc_fslmcflirt_outliers.csh

DESCRIPTION:
	This function 
	  1) uses mcflirt to do the motion correction. 
	  2) uses fsl_motion_outliers to obtain FDRMS and DVARS. 
	  3) uses DV_th and FD_th as threholds of DVARS and FDRMS to find outliers (high-motion frames). 
	     The default DV_th is 50; the default FD_th is 0.2. Frames either above FD_th or above DV_th 
	     will be regarded as outliers. 
	  4) 1 frame before, and 2 frames after each high-motion frame in (3) will also be labeled as 
	     outliers. The segments of low-motion frames lasting fewer than <discard_seg> are outliers as 
	     well. If -discard_seg not specified, the default is 5. 
	  5) discards the run which has more than <rm_run_th>% of frames of outliers.
	     In the output motion outlier text file, 0 means removing this frame, 1 means keeping this 
	     frame. 


REQUIRED ARGUMENTS:
	-s  <subject_id>           : subject's ID
	-d  <subject_dir>          : absolute path to <subject_id>. More specifically, processed data of 
	                             this subject will be in folder <subject_dir>/<subject_id>
	-bld  <bold_runs>          : bold run numbers, each number must have three digits. If this 
	                             subject has multiple runs, please use space as delimiter between 
	                             two run numbers (e.g. -bld "002 003"). NOTE: quote sign is necessary.
	-BOLD_stem  <BOLD_stem>    : specify the stem of input file. (e.g. if input file is 
	                             Sub0001_bld002_rest_skip4_stc.nii.gz then the stem of this file is 
	                             <BOLD_stem> = _rest_skip4_stc). This input file is assumed to be 
	                             stored in <subject_dir>/<subject_id>/bold/<run_number>.
	                             
OPTIONAL ARGUMENTS:
	-nframe  <nframe>          : perform mcflirt with the n-th frame. Default is <nframe> = 0, i.e. the 1st frame,
	-FD_th  <fd_th>            : FDRMS threshold for motion outlier detection, default <fd_th> = 0.2
	-DV_th  <dv_th>            : DVARS threshold for motion outlier detection, default <dv_th> = 50
	-force                     : force-update results. Existing results will be overwritten
	-discard-run  <rm_run_th>  : discard run which has more than <rm_run_th>% frames being outliers
	-rm-seg  <discard_seg>     : label the low-motion segments of data lasting fewer than <discard_seg> 
	                             contiguous frames as outliers.
	-spline_final              : interpolation method used in mcflirt, if the this option is used, the
								 interpolation method is spline, otherwise it is trilinear
	-help                      : help
	-version                   : version

OUTPUTS:
	(1) The NIFTI volume after motion correction:
	    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_mc.nii.gz
	    
	(2) A folder containing fsl_motion_outliers output files
		<subject_dir>/<subject_id>/bold/mc
		
	(3) mcflirt output files (motion parameters)
	    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_mc.par
	    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_mc_abs.rms
	    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_mc_abs_mean.rms
	    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_mc_rel.rms
	    <subject_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_mc_rel_mean.rms
	    
	    The last 4 files are made symbolic links to <subject_dir>/<subject_id>/qc folder.
	    
	(4) Motion outlier text file
	    In this file, each line is a number 0 or 1 corresponding to one frame. 0 represents this frame is 
	    censored (removed in analyses). 1 means this frame is uncensored (kept).
	    ${sub_dir}/${subject}/qc/${subject}_bld${run_folder}_FDRMS0.2_DVARS50_motion_outliers.txt
	
	(5) QC figures
	    
	    The traces of mcflirt estimated rotations (radians), translations (mm), and mean displacement (mm)
	    <subject_dir>/<subject_id>/qc/<subject_id>_bld<run_number><BOLD_stem>_mc_rot_trans_disp.png
	    
	    The scatter plot indicating the correlation between DVARS and FDRMS
	    <subject_dir>/<subject_id>/qc/<subject_id>_bld<run_number>_DVARS_FDRMS_correlation.png
	    
	    The traces of FDRMS and DVARS
	    <subject_dir>/<subject_id>/qc/<subject_id>_bld<run_number>_FDRMS<fd_th>_DVARS<dv_th>.png

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fslmcflirt_outliers.csh 
	-s Sub0033_Ses1 -d ~/storage/FMRI_preprocess -bld '002 003' -BOLD_stem _rest_skip4_stc -nframe 0 
	-FD_th 0.2 -DVARS 50 -discard-run 50 -rm-seg 5



