#! /bin/csh -f

#############################################
# Create masks and GLM regression
#############################################
# In this script, we: 
# 1) use CBIG_preproc_mask.csh to create masks (whole brain,wm,csf) 
# 2) use mcflirt results to create regressors (motion,whole brain,wm,csf)
# 3) use glm to regress out regressors for each run seperately
#
# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_regression.csh -s Su0033_Ses1 -d \
#	$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/100subjects_clustering/preproc_out -anat_s \
#	Sub0033_Ses1_FS -anat_d $CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/\
#   100subjects_clustering/recon_all -bld "002 003" -BOLD_stem _rest_skip4_stc_mc -REG_stem _rest_skip4_stc_mc_reg \
#   -MASK_stem _rest_skip4_stc_mc -erode_space anat -aCompCor -aCompCor_nPCs 5 -wm -wm_max_erode 3 -csf -csf_max_erode \
#	1 -motion12_itamar -per_run
#############################################
# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#############################################

set VERSION = '$Id: CBIG_preproc_regression.csh, v 1.0 2016/06/18 $'

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
set anat_s = ""
set anat_d = ""
set zpdbold = ""
set BOLD_stem = ""
set REG_stem = ""
set MASK_stem = ""
set OUTLIER_stem = ""

set whole_brain = 0		  # Default not create whole brain mask and regressor
set wm = 0				  # Default not create wm mask and regressor
set wm_max_erode = ""
set csf = 0				  # Default not create csf mask and regressor
set csf_max_erode = ""
set erode_space = "anat"  # Default space for wm and csf masks erosion is functional space
set motion12_itamar = 0	  # Default not create motion regressor
set aCompCor = 0          # Default not perform aCompCor
set nPCs = 5              # Default number of PCs of aCompCor is 5
set mt_diff_flag = 1;     # Default to include derivatives of motion parameters
set aCompCor_diff_flag = 0 # Default to not include derivatives of aCompCor PCs
set other_diff_flag = 1;   # Default to include derivatives of other regressors
set polynomial_fit = 1    # Default demean and remove the linear trend
set per_run = 0			  # Default regress out all runs jointly
set censor = 0			  # Default use all frames to do regression
set detrend_method = "detrend" # Default use detrend method when creating motion regressor
set force = 0

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

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
set LF = $subject_dir/$subject/logs/CBIG_preproc_regress.log
if( -e ${LF} ) then
    rm $LF
endif
touch $LF
echo "[Regression]: logfile = $LF"
echo "Regression" >> $LF
echo "[CMD]: CBIG_preproc_regress.csh $cmdline"   >>$LF

# boldfolder: directory of /bold
set boldfolder = "$subject_dir/$subject/bold"
echo "[Regression]: boldfolder = $boldfolder" |& tee -a $LF
echo "[Regression]: zpdbold = $zpdbold" |& tee -a $LF


#############################################
# Create masks (whole brain,wm,csf)
#############################################

echo "=======================Create masks (wb,wm,vent) for regressors=======================" |& tee -a $LF

if ( "$erode_space" == "func" ) then
	if ( $wm == 1 && "$wm_max_erode" == "" ) then
		set wm_max_erode = 1
		echo "WARNING: WM mask will be created but -wm_max_erode is not passed in. Maximal erosion " |& tee -a $LF
		echo "will be set to 1. (Erode in functional space.)" |& tee -a $LF
	endif
	
	if ( $csf == 1 && "$csf_max_erode" == "" ) then
		set csf_max_erode = 0
		echo "WARNING: Ventricles mask will be created but -csf_max_erode is not passed in. Maximal " |& tee -a $LF
		echo "erosion will be set to 0. (Erode in functional space.)" |& tee -a $LF
	endif
	
else
	if ( $wm == 1 && "$wm_max_erode" == "" ) then
		set wm_max_erode = 3
		echo "WARNING: WM mask will be created but -wm_max_erode is not passed in. Maximal erosion will " |& tee -a $LF
		echo "be set to 3. (Erode in anatomical space.)" |& tee -a $LF
	endif
	
	if ( $csf == 1 && "$csf_max_erode" == "" ) then
		set csf_max_erode = 1
		echo "WARNING: Ventricles mask will be created but -csf_max_erode is not passed in. Maximal " |& tee -a $LF
		echo "erosion will be set to 1. (Erode in anatomical space.)" |& tee -a $LF
	endif
endif

set cmd = "$root_dir/CBIG_preproc_create_mask.csh -s $subject -d $subject_dir -anat_s $anat -anat_d $anat_dir -bld "
set cmd = "$cmd '$zpdbold' -REG_stem $REG_stem -MASK_stem $MASK_stem"

if($whole_brain) then
	set cmd = "$cmd -whole_brain"
endif
if($wm) then
	set cmd = "$cmd -wm -wm_max_erode $wm_max_erode"
endif
if($csf) then
	set cmd = "$cmd -csf -csf_max_erode $csf_max_erode"
endif
if($wm || $csf || $aCompCor) then
	set cmd = "$cmd -erode_space $erode_space"
endif
if($aCompCor) then
	set cmd = "$cmd -aCompCor"
endif
if($force) then
	set cmd = "$cmd -force"
endif
echo $cmd |& tee -a $LF
eval $cmd 
echo "=======================Creating masks done!=======================" |& tee -a $LF



#############################################
# Create all regressors (mc,wb,wm,vent) 
#############################################

# check if matlab exists
set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
if ($status) then
	echo "ERROR: could not find matlab"
	exit 1;
endif

# create bold/regression folder
set regress_folder = "$boldfolder/regression"

if (! -e $regress_folder) then
    mkdir -p $regress_folder
endif

# create fMRI list, mc_par list, regressors list, censor list, resid list 
set fMRI_list = "$regress_folder/fMRI_list.txt" 
set mc_par_list = "$regress_folder/mc_par_list.txt"
set mc_regressor_list = "$regress_folder/mc_regressor_list.txt"
set ROI_regressors_list = "$regress_folder/ROI_regressors_list.txt"
set CompCor_regressors_list = "$regress_folder/aCompCor_regressors_list.txt";
set all_regressors_list = "$regress_folder/all_regressors_list.txt"
set censor_list = "$regress_folder/censor_list.txt"
set resid_list = "$regress_folder/resid_list.txt"

if (-e $fMRI_list) then
	rm $fMRI_list
endif
touch $fMRI_list

if (-e $mc_par_list) then
	rm $mc_par_list
endif
touch $mc_par_list

if (-e $mc_regressor_list) then
	rm $mc_regressor_list
endif
touch $mc_regressor_list

if (-e $ROI_regressors_list) then
	rm $ROI_regressors_list
endif
touch $ROI_regressors_list

if (-e $CompCor_regressors_list) then
	rm $CompCor_regressors_list
endif
touch $CompCor_regressors_list

if (-e $all_regressors_list) then
	rm $all_regressors_list
endif
touch $all_regressors_list

if (-e $censor_list) then
	rm $censor_list
endif
touch $censor_list

if (-e $resid_list) then
	rm $resid_list
endif
touch $resid_list


set mc_regressor_exist_flag = 1
set ROI_regressor_exist_flag = 1
set CompCor_regressor_exist_flag = 1
set resid_exist_flag = 1
set all_regressor_exist_flag = 1

# write file names to the lists and check whether output files already exist

# If $per_run is off, $CompCor_regressors_list should only contain one line 
# (only one CompCor regressors text file for merged runs)
if($per_run == 0) then
	set CompCor_regressor_file = "$regress_folder/${subject}_aCompCor_regressor.txt"
	
	if(! -e $CompCor_regressor_file) then
		set CompCor_regressor_exist_flag = 0
	endif
	
	echo $CompCor_regressor_file >> $CompCor_regressors_list
endif

# other regressors lists and final output
foreach curr_bold ($zpdbold)
	# file name of input fMRI volume
	set fMRI_file = "$boldfolder/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz"
	
	# file name of motion parameters from mcflirt
	set mc_par_file = `find $boldfolder/$curr_bold/*_mc.par`
	
	# file name of output motion regressors
	set mc_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}_mc_regressor.txt"
	
	# file name of output wb, wm, csf mean timeseries
	set ROI_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}"
	if ($whole_brain) then
		set ROI_regressor_file = "${ROI_regressor_file}_WB"
	endif
	if ($wm) then
		set ROI_regressor_file = "${ROI_regressor_file}_WM"
	endif
	if ($csf) then
		set ROI_regressor_file = "${ROI_regressor_file}_CSF"
	endif
	set ROI_regressor_file = "${ROI_regressor_file}_regressors.txt"
	
	# file name of output aCompCor regressors
	if ($per_run == 1) then
		set CompCor_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}_aCompCor_regressor.txt"
	endif
	
	set all_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}_all_regressors.txt"
	
	# file name of motion outliers
	if ( $censor == 1 ) then
		set censor_file = "$subject_dir/$subject/qc/$subject"_bld"${curr_bold}${OUTLIER_stem}"
	else
		set censor_file = ""
	endif
	
	# file name of output fMRI volume
	if ( $censor == 1 ) then
		set resid_file =  "$boldfolder/$curr_bold/$subject"_bld"${curr_bold}${BOLD_stem}_residc.nii.gz"
	else
		set resid_file =  "$boldfolder/$curr_bold/$subject"_bld"${curr_bold}${BOLD_stem}_resid.nii.gz"
	endif
	
	# check if regressors and output fMRI volume existance
	if (! -e $mc_regressor_file) then
		set mc_regressor_exist_flag = 0
	endif
	if (! -e $ROI_regressor_file) then
		set ROI_regressor_exist_flag = 0
	endif
	if ($per_run == 1) then
		if (! -e $CompCor_regressor_file) then
			set CompCor_regressor_exist_flag = 0
		endif
	endif
	if (! -e $all_regressor_file) then
		set all_regressor_exist_flag = 0
	endif
	if (! -e $resid_file) then
		set resid_exist_flag = 0
	endif
	
	# write file names into lists for each run
	echo $fMRI_file >> $fMRI_list
	echo $mc_par_file >> $mc_par_list
	echo $mc_regressor_file >> $mc_regressor_list
	echo $ROI_regressor_file >> $ROI_regressors_list
	if ($per_run == 1) then
		echo $CompCor_regressor_file >> $CompCor_regressors_list
	endif
	echo $all_regressor_file >> $all_regressors_list
	echo $censor_file >> $censor_list
	echo $resid_file >> $resid_list
end

# get mask file name
if ($whole_brain) then
	set wb_mask = "$boldfolder/mask/$subject.brainmask.bin.nii.gz" 
else
	set wb_mask = ""
endif
if ($wm) then
	set wm_mask = "$boldfolder/mask/$subject.func.wm.nii.gz"
else
	set wm_mask = ""
endif
if ($csf) then	
	set csf_mask = "$boldfolder/mask/$subject.func.ventricles.nii.gz"
else
	set csf_mask = ""
endif

# create mc regressor
echo "=======================Create mc regressor=======================" |& tee -a $LF
if ($motion12_itamar == 1) then
	if ($mc_regressor_exist_flag == 0) then
		set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r '"' 'addpath(fullfile('"'"$root_dir"'"\
','"'"utilities"'"'))'; CBIG_preproc_create_mc_regressors "'"$mc_par_list"'" "'"$mc_regressor_list"'" \
"'"$detrend_method"'" "'"$mt_diff_flag"'"; exit '"' );
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		
		foreach mc_regressor (${mc_regressor_list})
			if ( ! -e ${mc_regressor} ) then
				echo "ERROR: Motion regressor file ${mc_regressor} is not produced." |& tee -a $LF
				exit 1
			endif
		end
	else
		echo "mc regressor already exists." |& tee -a $LF
	endif
endif

# create ROI regressors
echo "=======================Create wb, wm, csf regressors=======================" |& tee -a $LF 
if ($whole_brain || $wm || $csf) then
	if ($ROI_regressor_exist_flag == 0) then
		set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r '"' 'addpath(fullfile('"'"$root_dir"'"\
','"'"utilities"'"'))'; CBIG_preproc_create_ROI_regressors "'"$fMRI_list"'" "'"$ROI_regressors_list"'" "'"$wb_mask"'" \
"'"$wm_mask"'" "'"$csf_mask"'" "'"$other_diff_flag"'"; exit '"' );
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		
		foreach ROI_regressor (${ROI_regressors_list})
			if ( ! -e ${ROI_regressor} ) then
				echo "ERROR: ROI regressor file ${ROI_regressor} is not produced." |& tee -a $LF
				exit 1
			endif
		end
	else
		echo "ROIs' mean regressors already exist." |& tee -a $LF
	endif
else
	echo "There are no mean wb, mean wm, or mean csf regressors."
endif

# create aCompCor regressors
echo "======================= Create aCompCor regressors =======================" |& tee -a $LF
if ( $aCompCor ) then
	if ($CompCor_regressor_exist_flag == 0) then
		set wm_csf_mask_comb = "$boldfolder/mask/${subject}.func.wm.vent.nii.gz"
		if(! -e $wm_csf_mask_comb) then
			echo "ERROR: The WB + ventricle mask does not exist." |& tee -a $LF
			exit 1;
		endif
		
		# The path needs to be changed after moved the matlab function in git repo
		set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r '"' 'addpath(fullfile('"'"$root_dir"'"','\
"'"utilities"'"'))'; CBIG_preproc_aCompCor_multipleruns "'"$fMRI_list"'" "'"$wm_csf_mask_comb"'" \
"'"$CompCor_regressors_list"'" "'"$nPCs"'" "'"$per_run"'" "'"$aCompCor_diff_flag"'" ; exit '"' );
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		
		foreach CompCor_regressor ($CompCor_regressors_list)
			if ( ! -e $CompCor_regressor ) then
				echo "ERROR: aCompCor regressor file ${CompCor_regressor} is not produced." |& tee -a $LF
				exit 1
			endif
		end
	else
		echo "aCompCor regressors already exist." |& tee -a $LF
	endif
else
	echo "Do not need to create aCompCor regressors." |& tee -a $LF
endif

# merge mc, ROI, aCompCor regressors
echo "=======================Merge all regressors=======================" |& tee -a $LF 
if (($all_regressor_exist_flag == 0) || ($force == 1)) then
	set i = 1
	foreach curr_bold ($zpdbold)
		set mc_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}_mc_regressor.txt"
		
		set ROI_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}"
		if ($whole_brain) then
			set ROI_regressor_file = "${ROI_regressor_file}_WB"
		endif
		if ($wm) then
			set ROI_regressor_file = "${ROI_regressor_file}_WM"
		endif
		if ($csf) then
			set ROI_regressor_file = "${ROI_regressor_file}_CSF"
		endif
		set ROI_regressor_file = "${ROI_regressor_file}_regressors.txt"
		
		set tmp1 = "$regress_folder/all_regressors_tmp.txt"
		set all_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}_all_regressors.txt"
		
		if ($motion12_itamar == 1 && ($wm == 1 || $csf == 1)) then
			paste -d " " $mc_regressor_file $ROI_regressor_file > $tmp1
		else if ($motion12_itamar == 1) then
			paste $mc_regressor_file > $tmp1
		else
			paste $ROI_regressor_file > $tmp1
		endif
		
		# pick aCompCor regressor
		if ($aCompCor) then
			if ($per_run == 1) then
				set CompCor_regressor_file = "$regress_folder/$subject"_bld"${curr_bold}_aCompCor_regressor.txt"
				paste -d " " $tmp1 $CompCor_regressor_file > $all_regressor_file
			else
				set CompCor_regressor_file = "$regress_folder/${subject}_aCompCor_regressor.txt"
				set tmp2 = "$regress_folder/aCompCor_regressor_tmp.txt"
				set func_name = "$boldfolder/$curr_bold/$subject"_bld"$curr_bold$BOLD_stem.nii.gz"
				set nframes = `fslval $func_name dim4`
				set start = `expr $i \* $nframes - $nframes + 1`
				set ending = `expr $i \* $nframes`
				sed -n "${start},${ending}p" $CompCor_regressor_file > $tmp2
				paste -d " " $tmp1 $tmp2 > $all_regressor_file
				
				#rm $tmp2
			endif
		else
			paste $tmp1 > $all_regressor_file
		endif
		
		#rm $tmp1
		set i = `expr $i + 1`
	end
else
	echo "All regressor already exists." |& tee -a $LF
endif

# if censor option is off
if ($censor == 0) then
	set censor_list = ""
endif

# start doing GLM regression
echo "=======================Use glm to regress out all regressors=======================" |& tee -a $LF
if ($motion12_itamar || $whole_brain || $wm || $csf || $aCompCor) then
	if (($resid_exist_flag == 0) || ($force == 1)) then
		set cmd = ( $MATLAB -nojvm -nodesktop -nodisplay -nosplash -r '"' 'addpath(fullfile('"'"$root_dir"'"','\
"'"utilities"'"'))'; CBIG_glm_regress_vol "'"$fMRI_list"'" "'"$resid_list"'" "'"$all_regressors_list"'" \
"'"$polynomial_fit"'" "'"$censor_list"'" "'"$per_run"'"; exit '"' )
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
	else
		echo "Resid already exists." |& tee -a $LF
	endif
else
	echo "No regressor, quit!" |& tee -a $LF
	exit 1;
endif
echo "=======================Regression done!=======================" |& tee -a $LF

#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "=======================Git: Last Commit of Current Function =======================" |& tee -a $LF
	pushd ${CBIG_CODE_DIR}
	git log -1 -- ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_regression.csh \
	>> $LF
	popd
endif

echo "*********************************************************************" |& tee -a $LF
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
			
		#anatomical name
		case "-anat_s":
			if ( $#argv == 0 ) goto arg1err;
			set anat = $argv[1]; shift;
			breaksw	
			
		#anatomical path
		case "-anat_d":
			if ( $#argv == 0 ) goto arg1err;
			set anat_dir = $argv[1]; shift;
			breaksw	
			
		#bold run number
		case "-bld"
			if ( $#argv == 0 ) goto arg1err;
			set zpdbold = ($argv[1]); shift;
			breaksw
			
		#input file stem
		case "-BOLD_stem"
			if ( $#argv == 0 ) goto arg1err;
			set BOLD_stem = $argv[1]; shift;
			breaksw
			
		#registration file stem
		case "-REG_stem"
			if ( $#argv == 0 ) goto arg1err;
			set REG_stem = $argv[1]; shift;
			breaksw
			
		#file stem used to create masks
		case "-MASK_stem"
			if ( $#argv == 0 ) goto arg1err;
			set MASK_stem = $argv[1]; shift;
			breaksw
			
		#file stem used to get censor frames
		case "-OUTLIER_stem"
			if ( $#argv == 0 ) goto arg1err;
			set OUTLIER_stem = $argv[1]; shift;
			breaksw

		#create whole brain mask and regressor 			
		case "-whole_brain":			
			set whole_brain = 1; 
			breaksw	
			
		#create wm mask and regressor
		case "-wm":			
			set wm = 1; 
			breaksw
		
		# maximal erode times for wm mask
		case "-wm_max_erode":
			if ( $#argv == 0 ) goto arg1err;
			set wm_max_erode = $argv[1]; shift;
			breaksw
			
		#create csf mask and regressor
		case "-csf":
			set csf = 1; 
			breaksw
			
		# maximal erode times for ventricles mask
		case "-csf_max_erode":
			if ( $#argv == 0 ) goto arg1err;
			set csf_max_erode = $argv[1]; shift;
			breaksw
		
		# erode space (func or anat):
		case "-erode_space":
			if ( $#argv == 0 ) goto arg1err;
			set erode_space = $argv[1]; shift;
			breaksw
			
		#create motion regressor
		case "-motion12_itamar":
			set motion12_itamar = 1; 
			breaksw
			
		# aCompCor
		case "-aCompCor":
			set aCompCor = 1;
			breaksw
			
		# number of PCs for aCompCor
		case "-aCompCor_nPCs":
			if ( $#argv == 0 ) goto arg1err;
			set nPCs = $argv[1]; shift;
			breaksw
			
		# whether derivatives of motion regressors will be included
		case "-no_mt_diff":
			set mt_diff_flag = 0;
			breaksw
			
		# whether derivatives of aCompCor regressors will be included
		case "-use_aCompCor_diff":
			set aCompCor_diff_flag = 1;
			breaksw
			
		# whether derivatives of other regressors will be included
		case "-no_other_diff":
			set other_diff_flag = 0;
			breaksw
			
		#add constant or linear detrend regressor
		case "-polynomial_fit"
			if ( $#argv == 0 ) goto arg1err;
			set polynomial_fit = $argv[1]; shift;
			breaksw

		#regress out the regressors for each run seperately	
		case "-per_run":
			set per_run = 1;
			breaksw
			
		#do regression without the censor frames, then apply the beta to all frames 	
		case "-censor"
			set censor = 1; 
			breaksw
			
		#use different detrend methods to detrend the motion regressors
		case "-detrend_method"
			if ( $#argv == 0 ) goto arg1err;
			set detrend_method = $argv[1]; shift;
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
if ( "$anat" == "" ) then
	echo "ERROR: anatomical folder not specified"
	exit 1;
endif
if ( "$anat_dir" == "" ) then
	echo "ERROR: anatomical directory not specified"
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
if ( "$REG_stem" == "" ) then
	echo "ERROR: registration file stem not specified"
	exit 1;
endif		
if ( "$MASK_stem" == "" ) then
	echo "ERROR: mask file stem not specified"
	exit 1;
endif

if ( "$erode_space" != "func" && "$erode_space" != "anat" ) then
	echo "ERROR: Wrong input for -erode_space. Please input func or anat."
	exit 1;
endif

if ( $censor == 1 && "$OUTLIER_stem" == "" ) then
	echo "ERROR: OUTLIER_stem not specified"
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
	CBIG_preproc_regression.csh
	
DESCRIPTION:
	Regress out regressors related to noise from fMRI volumes. To execute this function, we assume T1-T2* 
	registration of this subject has been finished.
	
	This function
		1) use CBIG_preproc_create_mask.csh to create masks (whole brain, wm, csf) 
		2) create regressors (motion, whole brain, wm, csf) 
		3) use glm to regress out regressors (if -per_run is used, regress each run separately; otherwise
		   regress all runs jointly).

REQUIRED ARGUMENTS:
	-s  <subject_id>                  : subject's ID
	-d  <subject_dir>                 : absolute path to <subject_id>. All preprocessed results of this subject 
	                                    will be at <subject_dir>/<subject_id>
	-anat_s  <anat_s>                 : anatomical folder for the subject (recon-all result)
	-anat_d  <anat_dir>               : absolute path to anatomical folder of this subject, i.e. recon-all 
	                                    results are in <anat_dir>/<anat_s>
	-bld  <bold_runs>                 : bold run numbers of this subject, each number should be three digits. 
	                                    If this subject has multiple runs, please use space as delimiter between 
	                                    two run numbers (e.g. -bld "002 003"). NOTE: quote sign is necessary.
	-BOLD_stem  <BOLD_stem>           : the stem of input file. E.g. if input file is 
	                                    Sub0001_bld002_rest_skip_stc_mc.nii.gz, then <BOLD_stem> = _rest_skip4_stc_mc. 
	                                    This input file is assumed to be in 
	                                    <subject_dir>/<subject_id>/bold/<run_number>.
	-REG_stem  <REG_stem>             : the stem of T1-T2* registration file. E.g. if the registration file 
	                                    is Sub0001_Ses1_bld002_rest_skip4_stc_mc_reg.dat, then 
	                                    <REG_stem> = _rest_skip4_stc_mc_reg. The registration file is assumed to be 
	                                    stored in <subject_dir>/<subject_id>/bold/<run_number>.
	                                    <subject_dir>/<subject_id>/bold/<run_number>.
	-MASK_stem  <MASK_stem>           : the stem of input file used to create the mask. E.g. if the template to 
	                                    create the masks is Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz, then 
	                                    <MASK_stem> = _rest_skip4_stc_mc. Notice that the choice of <MASK_stem> 
	                                    and <REG_stem> should be reasonable and consistent. For instance, the 
	                                    user cannot apply the registration file before motion correction on 
	                                    the template volume after motion correction.
	-OUTLIER_stem  <OUTLIER_stem>     : the stem of outlier text file. The number of lines in this file is equal to the 
	                                    total number of frames. Each line is a number of 0 or 1 (0-censored, 1-keep). 
	                                    If the outlier file is Sub0001_bld002_FDRMS0.2_DVARS50_motion_outliers.txt, then 
	                                    <OUTLIER_stem> is _FDRMS0.2_DVARS50_motion_outliers.txt. The outlier file is 
	                                    assumed to be stored in <subject_dir>/<subject_id>/qc.
	                                    
OPTIONAL ARGUMENTS:
	-whole_brain                      : if this option is used, create whole brain mask and regressor
	-wm                               : if this option is used, create white matter mask and regressor
	-wm_max_erode  <wm_max_erode>     : the maximal voxels to be eroded for WM mask. If erode in functional space, 
	                                    default is 1; if erode in anatomical space, default is 3.
	                                    All WM masks eroded from 0 to <wm_max_erode> will be created.
	                                    If erode in functional space, the final mask for WM is eroded once. 
	                                    If the WM mask and ventricles mask are eroded in anatomical space, the final 
	                                    mask for WM is the smallest but has at least 100 voxels.
	-csf                              : if this option is used, create ventricle mask and regressor
	-csf_max_erode  <csf_max_erode>   : the maximal voxels to be eroded for ventricles mask. If erode in functional 
	                                    space, default is 0; if erode in anatomical space, default is 1. 
	                                    All ventricles masks eroded from 0 to <csf_max_erode> will be created.
	                                    If erode in functional space, the final mask for ventricles is not eroded.
	                                    If erode in anatomical space, the final mask for ventricles is the one which 
	                                    is smallest but has at least 100 voxels.
	-erode_space  <erode_space>       : The space that the erosion of WM and ventricles masks will be done. 
	                                    Choose from "func" or "anat". Default is anat.
	-aCompCor                         : if this option is used, the WM mask and ventricles mask will be merged, and 
	                                    aCompCor regressors will be created.
	-aCompCor_nPCs  <nPCs>            : The number of PCs to be used as regressors for aCompCor. Default is 5.
	-motion12_itamar                  : if this option is used, create motion regressors. See 
	                                    CBIG_preproc_create_mc_regressors.m for more details.
	-no_mt_diff                       : If this option is used, derivatives of motion parameters will not be included. 
	                                    Default is to include the derivatives of motion parameters.
	-use_aCompCor_diff                : If this option is used, derivatives of aCompCor PCs will be included.
	                                    Default is to NOT include the derivatives of aCompCor PCs.
	-no_other_diff                    : If this option is used, the derivatives of other regressors except motion and  
	                                    aCompCor (GS, WM, CSF) will not be included. Default is to include the 
	                                    derivatives of other regressors.
	-polynomial_fit <polynomial_fit>  : -1/0/1 (default is 1). If polynomial_fit is set to -1, we add
	                                    nothing in regressor matrix. If polynomial_fit is set to 0, we prepend a
	                                    Mx1 vector [1,1,1...1]' to the regressor matrix and end up with a
	                                    Mx(K+1) matrix. If polynomial_fit is set to 1, a Mx1 vector
	                                    [1,1,1...1]' will be added into the first column and a Mx1 matrix
	                                    linspace(-1, 1, M)' will be added into the second column of the
	                                    regressor matrix. See CBIG_glm_regress_matrix.m for details.
	-per_run                          : If this option is used, regress out the regressors for each run separately. 
	                                    Default is to regress all runs jointly.
	-censor                           : fit beta coefficients without censored frames, then apply beta to all frames
	-force                            : update results, if output file exists then overwrite
	-help                             : help
	-version                          : version

OUTPUTS:
	(1) A folder <sub_dir>/<subject>/bold/regression
	    contains all regressor lists and regressors text files.
	    
	(2) The volume after regression
	    <sub_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_resid.nii.gz 
	    (if -censor is not used)
	    or
	    <sub_dir>/<subject_id>/bold/<run_number>/<subject_id>_bld<run_number><BOLD_stem>_residc.nii.gz 
	    (if -censor is used)

EXAMPLES:
	1.
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_regression.csh -s Sub0033_Ses1 \
	-d $CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/100subjects_clustering/preproc_out \
	-anat_s Sub0033_Ses1_FS -anat_d $CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/\
	100subjects_clustering/recon_all -bld "002 003" -BOLD_stem _rest_skip4_stc_mc -REG_stem _rest_skip4_stc_mc_reg \
	-MASK_stem _rest_skip4_stc_mc -whole_brain -wm -csf -motion12_itamar -per_run
	
	2.
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_regression.csh -s Sub0033_Ses1 
	-d $CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/100subjects_clustering/preproc_out \
	-anat_s Sub0033_Ses1_FS -anat_d $CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/\
	100subjects_clustering/recon_all -bld "002 003" -BOLD_stem _rest_skip4_stc_mc -REG_stem _rest_skip4_stc_mc_reg \
	-MASK_stem _rest_skip4_stc_mc -erode_space anat -aCompCor -aCompCor_nPCs 5 -wm -wm_max_erode 3 -csf -csf_max_erode \
	1 -motion12_itamar -per_run


