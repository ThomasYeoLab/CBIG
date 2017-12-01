#! /bin/csh -f
##########################################
# CBIG fMRI Preprocess
##########################################
# AUTHOR #################################
# RU(BY) KONG 
# 2016/06/09  
##########################################
##########################################
# In this script, we
# 1) Read in the configuration file to set the preprocess order
# 2) Obtain BOLD information from the fMRI nifti file list
# 3) Loop through the preprocess step
# 4) Preprocessed results will be saved in output directory
##########################################
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_fMRI_preprocess.csh, v 1.0 2016/06/09 $'


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

set subject = ""    # subject ID
set anat = ""       # recon-all folder
set anat_dir = ""   # path to recon-all folder
set output_dir = "" # output directory
set config = ""     # config file
set curr_stem = ""  # current stem
set zpdbold = ""    # bold runs (002 003)
set BOLD_stem = "_rest" # Initialize the BOLD_stem with _rest, the BOLD_stem is the input of each step and will be updated after each step
set REG_stem = ""   # stem of registration file produced by bbregister
set MASK_stem = ""  # stem of file used to create masks
set OUTLIER_stem="" # stem of file including censor frames
set nocleanup = 0   # default clean up intermediate files

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

##########################################
# Set preprocess log file and cleanup file
##########################################

mkdir -p $output_dir/$subject/logs
set LF = $output_dir/$subject/logs/CBIG_fMRI_preprocess.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
set cleanup_file = $output_dir/$subject/logs/cleanup.txt
if( -e $cleanup_file ) then
	rm $cleanup_file
endif
touch $cleanup_file
echo "**************************************************************************" >> $LF
echo "***************************CBIG fMRI Preprocess***************************" >> $LF
echo "**************************************************************************" >> $LF
echo "[LOG]: logfile = $LF" >> $LF
echo "[CMD]: CBIG_fMRI_preprocess.csh $cmdline"   >> $LF

##########################################
# Set env and git log file
##########################################

set TS = `date +"%Y-%m-%d_%H-%M-%S"`

# output env variables to env log
set env_log = $output_dir/$subject/logs/env.log
echo "$TS" >> $env_log
echo "***************************Env variable***************************" >> $env_log
env >> $env_log

# check if git exists
which git
if ($status) then
    echo "WARNING: could not find git, skip generating git log." >> $LF
else
	set git_log = $output_dir/$subject/logs/git.log
	echo "$TS" >> $git_log
	echo "***************************Git: Last Commit of Current Repo***************************" >> $git_log
	git -C ${CBIG_CODE_DIR} log -1 >> $git_log
endif

##########################################
# Read in Configuration file
##########################################

#filtering the comment line start with # and skip the blank lines
set config_clean = $output_dir/$subject/logs/CBIG_fMRI_preprocess.config
egrep -v '^#' $config | tr -s '\n' > $config_clean
set config = $config_clean

#print out the preprocessing order
echo "Verify your preprocess order:" >> $LF
foreach step ( "`cat $config`" )
	echo -n $step" => " >> $LF
end
echo "DONE!" >> $LF


##########################################
# Read in fMRI nifti file list
##########################################

#check fmri nifti file list columns
set lowest_numof_column = (`awk '{print NF}' $fmrinii_file | sort -nu | head -n 1`)
set highest_numof_column = (`awk '{print NF}' $fmrinii_file | sort -nu | tail -n 1`)
echo "lowest_numof_column = $lowest_numof_column" >> $LF
echo "highest_numof_column = $highest_numof_column" >> $LF
if ( $lowest_numof_column != 2 || $highest_numof_column != 2) then
	echo "[ERROR]: The input nifti file should only contain two columns!" >> $LF
	exit 1
endif
#check if there are repeating run numbers
set numof_runs_uniq = (`awk -F " " '{printf ("%03d\n", $1)}' $fmrinii_file | sort | uniq | wc -l`)
set zpdbold = (`awk -F " " '{printf ("%03d ", $1)}' $fmrinii_file`)
if ( $numof_runs_uniq != $#zpdbold ) then
	echo "[ERROR]: There are repeating bold run numbers!" >> $LF
	exit 1
endif

##########################################
# BOLD Information(read bold runs,output it into SUBJECT.bold file)
##########################################

echo "[BOLD INFO]: Number of runs: $#zpdbold" >> $LF 
echo "[BOLD INFO]: bold run $zpdbold" >> $LF
set boldname = (`awk -F " " '{printf($2" ")}' $fmrinii_file`)

#set output structure
@ k = 1
foreach curr_bold ($zpdbold)
	if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld$curr_bold$BOLD_stem.nii.gz" ) then
		mkdir -p $output_dir/$subject/bold/$curr_bold
		cp $boldname[$k] $output_dir/$subject/bold/$curr_bold/$subject"_bld$curr_bold$BOLD_stem.nii.gz"
	else
		echo "[BOLD INFO]: Input bold nifti file already exists !" >> $LF	
	endif
	echo "[BOLD INFO]: bold nifti file is $output_dir/$subject/bold/$curr_bold/$subject'_bld$curr_bold$BOLD_stem.nii.gz'" >> $LF
	@ k++
end
set Bold_file = $output_dir/$subject/logs/$subject.bold
if( -e $Bold_file ) then
	rm $Bold_file
endif
echo $zpdbold >> $Bold_file
echo "" >> $LF

##########################################
# Loop through each preprocess step
##########################################

foreach step ( "`cat $config`" )
	echo "" >> $LF
	
	#grep current preprocess step and input flag
	set curr_flag = ""
	set curr_step = (`echo $step | awk '{printf $1}'`)
	echo "[$curr_step]: Start..." >> $LF
	
	#get all arguments of current step
	set inputflag = (`echo $step | awk -F " " '{printf NF}'`)
	if ( $inputflag != 1) then
		set curr_flag = ( `echo $step | cut -f 1 -d " " --complement` )
		echo "[$curr_step] curr_flag = $curr_flag" >> $LF
	endif
	
	set zpdbold = `cat $Bold_file`
	echo "[$curr_step]: zpdbold = $zpdbold"	
	
	##########################################
	# Preprocess step: Skip first n frames 
	##########################################
	
	if ( "$curr_step" == "CBIG_preproc_skip" ) then
		
		set cmd = "$root_dir/CBIG_preproc_skip.csh -s $subject -d $output_dir -bld '$zpdbold' -BOLD_stem $BOLD_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null
		
		#update stem
		if ( $inputflag != 1 ) then
			set curr_stem = ("skip"`echo $curr_flag | awk -F "-skip" '{print $2}' | awk -F " " '{print $1}'`)
		else
			set curr_stem = "skip4"
		endif
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
			
		
		#check existence of output
		foreach curr_bold ($zpdbold)
		
		#put output into cleanup file
		echo $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz >> $cleanup_file
		
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_skip fail!" >> $LF
				exit 1
			endif
		end

	##########################################
	# Preprocess step: slice time correction 
	##########################################
		
	else if ( "$curr_step" == "CBIG_preproc_fslslicetimer" ) then
	    
	    set cmd = "$root_dir/CBIG_preproc_fslslicetimer.csh -s $subject -d $output_dir -bld '$zpdbold' -BOLD_stem $BOLD_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null 	
		
		#update stem	
		set curr_stem = "stc"
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
		
		#check existence of output
		foreach curr_bold ($zpdbold)
		
		#put output into cleanup file
		echo $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz >> $cleanup_file
		
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_fslslicetimer fail!" >> $LF
				exit 1
			endif
		end

	##########################################
	# Preprocess step: motion correction and detect outliers using FDRMS & DVARS as metric 
	##########################################		
		
	else if ( $curr_step == "CBIG_preproc_fslmcflirt_outliers" ) then
		
		set cmd = "$root_dir/CBIG_preproc_fslmcflirt_outliers.csh -s $subject -d $output_dir -bld '$zpdbold' -BOLD_stem $BOLD_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null
		
		#update stem 	
		set curr_stem = "mc"
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
		
		set FD_stem = (`echo $curr_flag | awk -F "-FD_th" '{print $2}' | awk -F " " '{print $1}'`)
		set DV_stem = (`echo $curr_flag | awk -F "-DV_th" '{print $2}' | awk -F " " '{print $1}'`)
		
		#Both FDRMS, DVARS threshold are given
		if ( "$FD_stem" != "" && "$DV_stem" != "") then 	
			set OUTLIER_stem = "_FDRMS${FD_stem}_DVARS${DV_stem}_motion_outliers.txt"
			echo "[$curr_step]: FDRMS threshold = $FD_stem"  >> $LF
			echo "[$curr_step]: DVARS threshold = $DV_stem"  >> $LF
		endif
		#Only FDRMS threshold given, use DVARS threshold as default: 50
		if ( "$FD_stem" != "" && "$DV_stem" == "") then 	
			set OUTLIER_stem = "_FDRMS${FD_stem}_DVARS50_motion_outliers.txt"
			echo "[$curr_step]: FDRMS threshold = $FD_stem"  >> $LF
			echo "[$curr_step]: DVARS threshold set as default: 50"  >> $LF
		endif
		#Only DVARS threshold given, use FDRMS threshold as default: 0.2
		if ( "$FD_stem" == "" && "$DV_stem" != "") then 	
			set OUTLIER_stem = "_FDRMS0.2_DVARS${DV_stem}_motion_outliers.txt"
			echo "[$curr_step]: FDRMS threshold set as default: 0.2"  >> $LF
			echo "[$curr_step]: DVARS threshold = $DV_stem"  >> $LF
		endif
		#FDRMS and DVARS threshold not given, use FDRMS threshold as default: 0.2, use DVARS threshold as default: 50
		if ( "$FD_stem" == "" && "$DV_stem" == "") then 	
			set OUTLIER_stem = "_FDRMS0.2_DVARS50_motion_outliers.txt"
			echo "[$curr_step]: FDRMS threshold set as default: 0.2"  >> $LF
			echo "[$curr_step]: DVARS threshold set as default: 50"  >> $LF
		endif
		echo "[$curr_step]: OUTLIER_stem = $OUTLIER_stem" >> $LF
		
        # if no run was left, then give a warning and exit the preprocessing 
	    set zpdbold = `cat $Bold_file`
        if ( "$zpdbold" == "" ) then
            echo "[WARNING]: There is no bold run left after discarding runs with high motion." >> $LF
            echo "Preprocessing Completed!" >> $LF
            exit 1
        endif
         
		#check existence of output
		foreach curr_bold ($zpdbold)
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_fslmcflirt_outliers fail!" >> $LF
				exit 1
			endif
		end
	
	##########################################
	# Preprocess step: function-anatomical registration
	##########################################
	
	else if ( "$curr_step" == "CBIG_preproc_bbregister" ) then
		
		set cmd = "$root_dir/CBIG_preproc_bbregister.csh -s $subject -d $output_dir -anat_s $anat -bld '$zpdbold' -BOLD_stem $BOLD_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null	
		
		#update stem
		set curr_stem = "reg"
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set REG_stem = $BOLD_stem
		set MASK_stem = $BOLD_stem
		set REG_stem = $REG_stem"_$curr_stem"
		
		#check existence of output
		foreach curr_bold ($zpdbold)
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$REG_stem.dat ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$REG_stem.dat can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_bbregister!" >> $LF
				exit 1
			endif
		end		

    ##########################################
	# Preprocess step: despiking 
	##########################################
		
	else if ( "$curr_step" == "CBIG_preproc_despiking" ) then
	    
	    set cmd = "$root_dir/CBIG_preproc_despiking.csh -s $subject -d $output_dir -bld '$zpdbold' -BOLD_stem $BOLD_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null 	
		
		#update stem	
		set curr_stem = "dspk"
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
		
		#check existence of output
		foreach curr_bold ($zpdbold)	
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_despiking fail!" >> $LF
				exit 1
			endif
		end	
		

	##########################################
	# Preprocess step: Motion Scrubbing and Interpolation
	##########################################
	
	else if ( $curr_step == "CBIG_preproc_censor" ) then
		
		# usage of -nocleanup option of censoring interpolation step is allowing the wrapper function
		if ( $nocleanup != 1) then            # needs to cleanup
			# if curr_flag contains "-nocleanup", we remove this option
			set curr_flag = `echo $curr_flag | sed 's/-nocleanup//'`
			echo "[$curr_step]: -nocleanup is not passed in wrapper function CBIG_fMRI_preprocess.csh. The intermediate censoring interpolation volume will be removed." >> $LF
		else                                  # does not need to cleanup
			# check if curr_flag contains -nocleanup
			set flag_ind = `echo $curr_flag | awk '{print index($0, "-nocleanup")}'`
			if ( $flag_ind == 0 ) then        # 0 means curr_flag does not contain "-nocleanup", add "-nocleanup"
				set curr_flag = "$curr_flag -nocleanup"
			endif
			echo "[$curr_step]: -nocleanup is passed in wrapper function CBIG_fMRI_preprocess.csh. The intermediate censoring interpolation volume will not be removed." >> $LF
		endif
		
		set cmd = "$root_dir/CBIG_preproc_censor.csh -s $subject -d $output_dir -anat_s $anat -anat_d $SUBJECTS_DIR -bld '$zpdbold' -BOLD_stem $BOLD_stem -REG_stem $REG_stem -OUTLIER_stem $OUTLIER_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null 
		
		#update stem
		set FD_th = (`echo $OUTLIER_stem | awk -F "FDRMS" '{print $2}' | awk -F "_" '{print $1}'`)
		set DV_th = (`echo $OUTLIER_stem | awk -F "DVARS" '{print $2}' | awk -F "_" '{print $1}'`)	
		set low_f = ( `echo $curr_flag | awk -F "-low_f" '{print $2}' | awk -F " " '{print $1}'` )
		set high_f = ( `echo $curr_flag | awk -F "-high_f" '{print $2}' | awk -F " " '{print $1}'` )
		set curr_stem = "interp_FDRMS${FD_th}_DVARS${DV_th}"
		if ( "$low_f" != "" && "$high_f" != "") then 
			set curr_stem = ${curr_stem}"_bp_"$low_f"_"$high_f
		endif
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
		 
		#check existence of output
		set zpdbold = `cat $Bold_file`
           
		foreach curr_bold ($zpdbold)
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_censor fail!" >> $LF
				exit 1
			endif
		end

	##########################################
	# Preprocess step: Bandpass/Lowpass/Highpass Filtering
	##########################################
	
	else if ( "$curr_step" == "CBIG_preproc_bandpass" ) then
		
		set cmd = "$root_dir/CBIG_preproc_bandpass_fft.csh -s $subject -d $output_dir -bld '$zpdbold' -BOLD_stem $BOLD_stem -OUTLIER_stem $OUTLIER_stem $curr_flag" 
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null 	
		
		#update stem
		set low_f = ( `echo $curr_flag | awk -F "-low_f" '{print $2}' | awk -F " " '{print $1}'` )
		set high_f = ( `echo $curr_flag | awk -F "-high_f" '{print $2}' | awk -F " " '{print $1}'` )
		#bandpass
		if ( "$low_f" != "" && "$high_f" != "") then 
			set curr_stem = "bp_"$low_f"_"$high_f
			echo "[$curr_step]: bandpass filtering.." >> $LF
			echo "[$curr_step]: low_f = $low_f"  >> $LF
			echo "[$curr_step]: high_f = $high_f"  >> $LF
		else
			echo "Please specify both low frequency and high frequency, e.g. -low_f 0.001 -high_f 0.08"
			exit 1
		endif
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
		
		#check existence of output
		foreach curr_bold ($zpdbold)
		
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_bandpass fail!" >> $LF
				exit 1
			endif
		end

	##########################################
	# Preprocess step: Regression ( Use ouput of last step as the MASK input )
	##########################################
	
	else if ( $curr_step == "CBIG_preproc_regress" ) then
		
		set cmd = "$root_dir/CBIG_preproc_regression.csh -s $subject -d $output_dir -anat_s $anat -anat_d $SUBJECTS_DIR -bld '$zpdbold' -BOLD_stem $BOLD_stem -REG_stem $REG_stem -MASK_stem $BOLD_stem -OUTLIER_stem $OUTLIER_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null 
		
		#update stem
		set censor_flag = ( `echo $curr_flag | grep "-censor"` )
		if ( $censor_flag == 1 ) then
			set curr_stem = "residc"
		else	
			set curr_stem = "resid"
		endif
		
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set BOLD_stem = $BOLD_stem"_$curr_stem"
		#check existence of output
		foreach curr_bold ($zpdbold)
		
		#put output into cleanup file
		echo $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz >> $cleanup_file
		
			if ( ! -e $output_dir/$subject/bold/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/bold/$curr_bold/${subject}_bld$curr_bold$BOLD_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_regress fail!" >> $LF
				exit 1
			endif
		end	

	##########################################
	# Preprocess step: Porjection to fsaverage surface (project to high resolution => smooth => downsample to low resolution)
	##########################################
	
	else if ( $curr_step == "CBIG_preproc_native2fsaverage" ) then
		
		set cmd = "$root_dir/CBIG_preproc_native2fsaverage.csh -s $subject -d $output_dir -anat_s $anat -anat_d $SUBJECTS_DIR -bld '$zpdbold' -BOLD_stem $BOLD_stem -REG_stem $REG_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null
		
		#update stem
		if ( $inputflag != 1 ) then
			set proj_mesh = ( `echo $curr_flag | awk -F "-proj" '{print $2}' | awk -F " " '{print $1}'` )
			set sm = ( `echo $curr_flag | awk -F "-sm" '{print $2}' | awk -F " " '{print $1}'` )
			set down_mesh = ( `echo $curr_flag | awk -F "-down" '{print $2}' | awk -F " " '{print $1}'` )
			set proj_res = `echo -n $proj_mesh | tail -c -1`
			if($proj_res == "e") then
				set proj_res = 7
			endif

			set down_res = `echo -n $down_mesh | tail -c -1`
			if($down_res == "e") then
				set down_res = 7;
			endif
		else
			set curr_stem = "fs6_sm6_fs5"
		endif
				
		set curr_stem = fs${proj_res}_sm${sm}_fs${down_res}
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set SURF_stem = $BOLD_stem"_$curr_stem"
		set FC_SURF_stem = ${BOLD_stem}_fs${proj_res}_sm${sm}
		
		#check existence of output
		foreach curr_bold ($zpdbold)
			if ( ! -e $output_dir/$subject/surf/lh.$subject"_bld"$curr_bold$SURF_stem.nii.gz || ! -e $output_dir/$subject/surf/rh.$subject"_bld"$curr_bold$SURF_stem.nii.gz) then
				echo "[ERROR]: file $output_dir/$subject/surf/${subject}_bld$curr_bold$SURF_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_native2fsaverage fail!" >> $LF
				exit 1
			endif
		end	
	
	##########################################
	# Preprocess step: Compute FC (functional connectivity) metrics
	##########################################
	else if ( $curr_step == "CBIG_preproc_FC_metrics" ) then
		
		# usage of -nocleanup option of censoring interpolation step is allowing the wrapper function
		if ( $nocleanup != 1) then            # needs to cleanup
			# if curr_flag contains "-nocleanup", we remove this option
			set curr_flag = `echo $curr_flag | sed 's/-nocleanup//'`
			echo "[$curr_step]: -nocleanup is not passed in wrapper function CBIG_fMRI_preprocess.csh. \
The intermediate files will be removed." >> $LF
		else                                  # does not need to cleanup
			# check if curr_flag contains -nocleanup
			set flag_ind = `echo $curr_flag | awk '{print index($0, "-nocleanup")}'`
			if ( $flag_ind == 0 ) then        # 0 means curr_flag does not contain "-nocleanup", add "-nocleanup"
				set curr_flag = "$curr_flag -nocleanup"
			endif
			echo "[$curr_step]: -nocleanup is passed in wrapper function CBIG_fMRI_preprocess.csh. \
			                    The intermediate files will not be removed." >> $LF
		endif
		
		set cmd = "$root_dir/CBIG_preproc_FCmetrics_wrapper.csh -s $subject -d $output_dir -bld '$zpdbold' "
		set cmd = "$cmd -BOLD_stem $BOLD_stem -SURF_stem $FC_SURF_stem -OUTLIER_stem $OUTLIER_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >& /dev/null
		
		# update stem & check existance of output
		if ( "$curr_flag" =~ *"-Pearson_r"* ) then
			set FC_metrics_stem = "${FC_SURF_stem}_all2all"
			if ( ! -e $output_dir/$subject/FC_metrics/Pearson_r/$subject$FC_metrics_stem.mat ) then
				echo "[ERROR]: file $output_dir/$subject/FC_metrics/Pearson_r/$subject$FC_metrics_stem.mat can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_FC_metrics fail!" >> $LF
				exit 1
			endif
		endif

	##########################################
	# Preprocess step: Porjection to MNI volume space (Project to FS1mm => MNI1mm => MNI2mm => Smooth)
	##########################################
	
	else if ( $curr_step == "CBIG_preproc_native2mni" ) then
		
		set cmd = "$root_dir/CBIG_preproc_native2mni.csh -s $subject -d $output_dir -anat_s $anat -anat_d $SUBJECTS_DIR -bld '$zpdbold' -BOLD_stem $BOLD_stem -REG_stem $REG_stem $curr_flag"
		echo "[$curr_step]: $cmd" >> $LF
		eval $cmd >&  /dev/null 	
		
		#update stem
		if ( $inputflag != 1 ) then
			set sm = ( `echo $curr_flag | awk -F "-sm" '{print $2}' | awk -F " " '{print $1}'` )
			set curr_stem = "FS1mm_MNI1mm_MNI2mm_sm"$sm
		else
			set curr_stem = "FS1mm_MNI1mm_MNI2mm_sm6"
		endif
		set final_mask = (`echo $curr_flag | awk -F "-final_mask" '{print $2}' | awk -F " " '{print $1}'`)
		if ( $final_mask != "" ) then
			set curr_stem = ${curr_stem}_finalmask
		endif
		echo "[$curr_step]: bold_stem = $curr_stem" >> $LF
		set VOL_stem = $BOLD_stem"_$curr_stem"
		#check existence of output
		foreach curr_bold ($zpdbold)
			if ( ! -e $output_dir/$subject/vol/$subject"_bld"$curr_bold$VOL_stem.nii.gz ) then
				echo "[ERROR]: file $output_dir/$subject/vol/${subject}_bld$curr_bold$VOL_stem.nii.gz can not be found" >> $LF
				echo "[ERROR]: CBIG_preproc_native2mni fail!" >> $LF
				exit 1
			endif
		end
		
	else
	
	##########################################
	# Preprocess step can not be recognized 
	##########################################	
	
		echo "ERROR: $curr_step can not be identified in our preprocessing step" >> $LF
		exit 1
		
	endif
			
end

#########################################
# echo successful message
#########################################
echo "Preprocessing Completed!" >> $LF

	
##########################################
# clean up intermediate files 
##########################################	
if ( $nocleanup != 1) then
foreach file (`cat $cleanup_file`)
	rm $file
end
endif
exit 1

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
		case "-fmrinii":
			if ( $#argv == 0 ) goto arg1err;
			set fmrinii_file = $argv[1]; shift;
			breaksw
			
		#anatomical ID
		case "-anat_s":
			if ($#argv == 0) goto arg1err;
			set anat = $argv[1]; shift;
			breaksw
			
		#directory to recon-all folder
		case "-anat_d":
			if ($#argv == 0) goto arg1err;
			set anat_dir = $argv[1]; shift;
			setenv SUBJECTS_DIR $anat_dir
			breaksw
			
		#output directory to save out preprocess results
		case "-output_d":
			if ( $#argv == 0 ) goto arg1err;
			set output_dir = $argv[1]; shift;
			breaksw	
				
		#configuration file		
		case "-config":
			if ( $#argv == 0 ) goto arg1err;
			set config = $argv[1]; shift;
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

##########################################
# Check Parameters
##########################################

check_params:
if ( "$subject" == "" ) then
	echo "ERROR: subject not specified"
	exit 1;
endif

if ( "$fmrinii_file" == "" ) then
	echo "ERROR: subject's fmri nifti file list not specified"
	exit 1;
endif

if ( "$anat" == "" ) then
	echo "ERROR: subject's anatomical ID not specified"
	exit 1;
endif

if ( "$anat_dir" != "$SUBJECTS_DIR" ) then
	echo "ERROR: subject's anatomical data directory doesn't match enviromental variable SUBJECTS_DIR!"
	echo "ERROR: anat_dir = $anat_dir"
	echo "ERROR: SUBJECTS_DIR = $SUBJECTS_DIR"
	exit 1;
endif

if ( "$output_dir" == "" ) then
	echo "ERROR: preprocess output directory not specified"
	exit 1;
endif

if ( "$config" == "" ) then
	echo "ERROR: subject's configuration file not specified"
	exit 1;
endif				
goto check_params_return;

##########################################
# ERROR message
##########################################
	
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
  
#####################################
# Help
#####################################
BEGINHELP

NAME:
	CBIG_fMRI_preprocess.csh

DESCRIPTION:
	The pipeline processes fMRI data and projects the data to 
	(1) FreeSurfer fsaverage5, fsaverage6 space
	(2) FreeSurfer nonlinear volumetric space
	(3) FSL MNI152 1mm, MNI152 2mm space.
        
	The pipeline proceeds sequentially as follows (default order), you can change the order and parameters 
	by changing the config file:
	(1) [CBIG_preproc_skip -skip 4] 
	    skips first 4 frames of resting data. 
	(2) [CBIG_preproc_fslslicetimer -slice_order <so_file>] 
	    does slice time correction using FSL slicetimer. If the user does not pass in the slice acquisition direction 
	    -direction <direction>, this step will use "Siemens" acquisition direction Superior-Inferior as default. 
	    If the direction is Right-Left, <direction> should be 1 representing x axis.
	    If the direction is Anterior-Posterior, <direction> should be 2 representing y axis.
	    If the direction is Superior-Inferior, <direction> should be 3 representing z axis.
	    If the user does not pass in the slice order file <so_file>, this step will use "Siemens" ordering as default. 
	    If the number of slices is odd, the ordering is 1, 3, 5, ..., 2, 4, 6, ...; 
	    if the number of slices is even, the ordering is 2, 4, 6, ..., 1, 3, 5, ....
	(3) [CBIG_preproc_fslmcflirt_outliers -FD_th 0.2 -DV_th 50 -discard-run 50 -rm-seg 5 -spline_final] 
	    does motion correction with spline interpolation and calculates Framewise Displacement and DVARS, 
	    then generates a vector indicating censored frames (1:keep 0:censored). This step throws away the
	    runs where the number of outliers are more than the threshold set by -discard-run option.
	(4) [CBIG_preproc_bbregister -intrasub_best] 
	    does bbregister with fsl initialization for each run and chooses the best run to generate registration file.
	(5) [CBIG_preproc_regress -whole_brain -wm -csf -motion12_itamar -detrend_method detrend -per_run -censor \
	    -polynomial_fit 1] 
	    regresses out motion, whole brain, white matter, ventricle, linear regressors for each run seperately. 
	    If the data have censored frames, this function will first estimate the beta coefficients ignoring the 
	    censored frames and then apply the beta coefficients to all frames to regress out those regressors.  
	(6) [CBIG_preproc_censor -nocleanup] 
	    removes (ax+b) trend of censored frames, then does censoring with interpolation. For interpolation method, 
	    please refer to (Power et al. 2014).
	(7) [CBIG_preproc_bandpass -low_f 0.009 -high_f 0.08 -detrend] 
	    does bandpass filtering with passband = [0.009, 0.08] (boundaries are included). This step applies FFT 
	    on timeseries and cuts off the frequencies in stopbands (rectanguluar filter), then performs inverse FFT 
	    to get the result.
	(8) [CBIG_preproc_native2fsaverage -proj fsaverage6 -down fsaverage5 -sm 6] 
	    projects fMRI to fsaverage6, smooths it with fwhm = 6mm and downsamples it to fsaverage5.
	(9) [CBIG_preproc_FC_metrics -Pearson_r -censor -lh_cortical_ROIs_file <lh_cortical_ROIs_file> -rh_cortical_ROIS_file \
	    <rh_cortical_ROIs_file>]
	    computes FC (functional connectivity) metrics based on both cortical and subcortical ROIs. The cortical ROIs 
	    can be passed in by -lh_cortical_ROIs and -rh_cortical_ROIs. The subcortical ROIs are 19 labels extracted 
	    from aseg in subject-specific functional space. This function will support for multiple types of FC metrics
	    in the future, but currently we only support static Pearson's correlation by using "-Pearson_r" flag. 
	    If "-censor" flag is used, the censored frames are ignored when computing FC metrics.
	(10) [CBIG_preproc_native2mni -down FSL_MNI_2mm -sm 6 -sm_mask <sm_mask> -final_mask <final_mask>] 
	    first, projects fMRI to FreeSurfer nonlinear space; second, projects from FreeSurfer nolinear space to 
	    FSL MNI 1mm space; third, downsamples from FSL MNI 1mm space to FSL MNI 2mm space; fourth, smooths it by 
	    fwhm = 6mm within <sm_mask>; and last, masks the result by <final_mask> to save disk space.
	
	There is a despiking function which is not used in the default config file:
	[CBIG_preproc_despiking]
	It uses AFNI 3dDespike to conduct despiking. This function can be used to replace censoring interpolation step, depended  
	on the requirement of users.

	Note: this pipeline assumes the user has finished FreeSurfer recon-all T1 preprocessing.
	
	This pipeline also assumes that no reorientation is needed. To decide whether reorientation is required, 
	load volumes on freeview. If the following statements are true, you do not need to reorient your data:
	- Scroll down coronal, 2nd voxel coordinate decreases
	- Scroll down sagittal, 1st voxel coordinate increases
	- Scroll down axial, 3rd voxel coordinate decreases
	If reorientation is needed, please refer to fslreorient2std command.
	 
	Please be aware that T1-T2* registration step must be done before CBIG_preproc_regress and CBIG_preproc_censor. 
	The latter two steps need registration information to create masks.
	
	To know how to do QC checks, please refer to README.md in the same folder as this script.

   
REQUIRED ARGUMENTS:
	-s  <subject>              : subject ID
	-fmrinii  <fmrinii>        : fmrinii text file including 2 columns, the 1st column contains all run numbers, 
	                            the 2nd column specify the absolute path to raw functional nifti files for corresponding run.
	                            An example file is here: 
	                            ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/example_fmrinii.txt
	                            Example of <fmrinii> content:
	                            002 /data/../Sub0015_bld002_rest.nii.gz
	                            003 /data/../Sub0015_bld003_rest.nii.gz

	-anat_s  <anat>            : FreeSurfer recon-all folder name of this subject (relative path)
	-anat_d  <anat_dir>        : specify anat directory to recon-all folder (full path), i.e. <anat_dir> contains <anat>
	-output_d  <output_dir>    : output directory to save out preprocess results (full path). This pipeline will create a 
	                             folder named <subject> under <output_dir>. All preprocessing results of this subject are 
	                             stored in <output_dir>/<subject>.
	-config  <config>          : configuration file
	                            An example file is here: 
	                            ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/example_config.txt
	                            Example of <config> content (Remind: this is not a full config file):
	                            ###CBIG fMRI preprocessing configuration file
	                            ###The order of preprocess steps is listed below
	                            CBIG_preproc_skip -skip 4
	                            CBIG_preproc_fslslicetimer -slice_order ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/example_slice_order.txt
	                            CBIG_preproc_fslmcflirt_outliers -FD_th 0.2 -DV_th 50 -discard-run 50 -rm-seg 5

	                            The symbol # in the config file also means comment, you can write anything you want if you 
	                            begin a line with #. Each line of the config file representing a function or step of our 
	                            preprocessing pipeline, the order of the step representing our preprocessing order, so it is 
	                            easy to change the order of your preprocessing according to changing the config file.
	                            In this config file, you can also specify the option of each function. For example, if you 
	                            want to skip first 4 frames of the fMRI data, you can add the option (-skip 4) behind the 
	                            CBIG_preproc_skip. For further details about these options, you can use option (-help) for 
	                            each function, such as (CBIG_preproc_skip -help).

OPTIONAL ARGUMENTS:
	-help                      : help
	-version                   : version
	-nocleanup                 : do not delete intermediate volumes

OUTPUTS: 
	CBIG_fMRI_preprocess.csh will create the directory <output_dir>/<subject> as specified in the options. Within the 
	<output_dir>/<subject> folder, there are multiple folders:

	1. surf folder contains the intermediate and final preprocessed fMRI data on the surface. 
		For example, 
		surf/lh.Sub0033_Ses1_bld002_rest_skip4_stc_mc_resid_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz 
		is bold data from run 002 ("bld002") of subject "Sub0033_Ses1" that has been projected to the left hemisphere ("lh"). 
		The remaining descriptors in the filename describes the ordering of the processing that has occurred. In particular,
		"rest" = resting state fmri
		"skip" = first four frames have been removed for T1 equilibrium
		"stc" = slice time correction
		"mc" = motion correction
		"resid" = regression of motion, whole brain, ventricular, white matter signals (standard fcMRI preprocessing)
		"interp_FDRMS0.2_DVARS50" = do interpolation for the censored frames defined by Framewise Displacement > 0.2,
		                            DVARS > 50, 
		"bp_0.009_0.08" = bandpass filtering with passband = [0.009, 0.08] (boundary inclusive).
		"fsaverage6" = data projected to fsaverage6 surface
		"sm6" = data smoothed with a FWHM = 6mm kernel on the surface
		"fsaverage5" = data downsampled to fsaverage5 surface

	2. vol folder contains the intermediate and final preprocessed fMRI data in the MNI152 and freesurfer nonlinear 
	   volumetric space.
		For example, 
		vol/Sub0033_Ses1_bld002_rest_skip4_stc_mc_resid_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_FS1mm_MNI1mm_MNI2mm_sm6.nii.gz 
		is bold data from run 002 ("bld002") subject "Sub0033_Ses1". The remaining descriptors in the filename describes the 
		ordering of the processing that has occurred. In particular,
		"rest" = resting state fmri
		"skip" = first four frames have been removed for T1 equilibrium
		"stc" = slice time correction
		"mc" = motion correction
		"resid" = regression of motion, whole brain, ventricular, white matter signals (standard fcMRI preprocessing)
		"interp_FDRMS0.2_DVARS50" = do interpolation for the censored frames defined by Framewise Displacement > 0.2, 
		                            DVARS > 50, 
		"bp_0.009_0.08" = bandpass filtering with passband = [0.009, 0.08] (boundary inclusive).
		"FS1mm" = projection of data to freesurfer nonlinear 1mm volumetric space
		"MNI1mm" = projection of data to MNI152 nonlinear 1mm volumetric space
		"MNI2mm" = downsampling of data to MNI152 nonlinear 2mm volumetric space
		"sm6" = data smoothed with a FWHM = 6mm kernel
		
	3. logs folder contains all log files for our preprocessing.
		CBIG_fMRI_preprocess.log contains the log info of CBIG_fMRI_preprocess.csh function, which is a wrapper script.
		Similarly, the name of the log file indicates the function, for example, CBIG_preproc_regress.log corresponds to the 
		function CBIG_preproc_regression.csh. Other log files: env.log includes all environment variables; git.log includes 
		the last git commit info; Sub0033_Ses1.bold contains the run numbers of this subject after censoring; cleanup.txt 
		includes all intermediate files that have been deleted, the user can use -nocleanup option to keep these volumes.
	   
	4. bold folder contains the intermediate files for each step.
		bold/002 folder contains all intermediate bold volumes of run 002.
		For example, Sub0033_Ses1_bld002_rest_skip4_stc_mc.nii.gz is the volume after skip -> slice-timing correction -> 
		motion correction

		bold/mask folder contains all the fMRI masks.
		For example, Sub0033_Ses1.func.ventricles.nii.gz means that it's a functional ventricle mask for the subject 
		Sub0033_Ses1; Sub0033_Ses1.brainmask.bin.nii.gz means that it's a binarized brainmask for subject Sub0033_Ses1.

		bold/regression folder contains all regressors and lists for the glm regression.
		For example, Sub0033_Ses1_bld002_all_regressors.txt means all regressors of subject Sub0033_Ses1, run 002.

		bold/mc folder contains the output files of fsl_motion_outliers, and some intermediate files when detecting 
		high-motion outliers. For example, Sub0033_Ses1_bld002_rest_skip4_stc_motion_outliers_DVARS is the text file of 
		DVARS value of each frame of Sub0033_Ses1, run 002;	Sub0033_Ses1_bld002_rest_skip4_stc_motion_outliers_FDRMS is the 
		text file of FDRMS value of each frame of Sub0033_Ses1, run 002;

	5. qc folder contains all the files that are useful for quality control.
		For example:
		CBIG_preproc_bbregister_intra_sub_reg.cost contains the number of bbregister cost in T1-T2* registration.
		Sub0033_Ses1_bld002_mc_abs.rms, Sub0033_Ses1_bld002_mc_abs_mean.rms, Sub0033_Ses1_bld002_mc_rel.rms, and 
		Sub0033_Ses1_bld002_mc_rel_mean.rms are motion parameters.
		Sub0033_Ses1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt contains the outlier labels of frames (1-keep, 0-censored)
		For introduction of more qc files, please refer to README.md in the same folder of this script.
		
	6. FC_metrics folder contains all files related to this subject's FC (functional connectivity) metrics.
	   It contains three subfolders currently"
	   FC_metrics/ROIs contains the 19 subcortical ROIs file;
	   FC_metrics/lists contains the input lists for corresponding matlab function;
	   FC_metrics/Pearson_r contains the static Pearson's correlation of this subject.
  
EXAMPLE:
	CBIG_fMRI_preprocess.csh -s Sub0033_Ses1 -output_d ./ -anat_s Sub0033_Ses1_FS -anat_d 
	/share/users/imganalysis/yeolab/data/GSP_release -fmrinii ./Sub0033_Ses1.fmrinii -config ./prepro.config
	
Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

