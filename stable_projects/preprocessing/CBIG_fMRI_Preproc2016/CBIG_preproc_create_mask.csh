#! /bin/csh -f

#############################################
# Create fMRI masks (wb,wm,gm,csf)
#############################################
# AUTHOR ####################################
# Nanbo Sun
# 2016/06/18  
#############################################
#############################################
# In this script, we: 
# 1) Create fMRI masks (wb,wm,gm,csf)
# Example: 
#	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_mask.csh -s Sub0001_Ses1 
#	-d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS 	-anat_d ~/storage/sMRI_preprocess -bld "002 003" 
#	-REG_stem _rest_skip4_stc_mc_reg -MASK_stem _rest_skip4_stc_mc -whole_brain -wm -csf -gm
#############################################
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_create_mask.csh, v 1.0 2016/06/18 $'

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
set anat = ""
set anat_dir = ""
set REG_stem = ""
set MASK_stem = ""
set whole_brain = 0       # Default not create whole_brain mask
set loose_whole_brain = 0 # Default not create loose_whole_brain mask
set wm = 0 			      # Default not create wm mask
set wm_max_erode = ""
set csf = 0			      # Default not create csf mask
set csf_max_erode = ""
set erode_space = "anat"  # Default space for wm and csf masks erosion is anatomical space
set gm = 0			      # Default not create gm mask
set aCompCor = 0
set force = 0		      # Default if file exist, skip the step

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

if ( $aCompCor == 1 ) then
	set wm = 1
	set csf = 1
endif

cd $subject_dir/$subject

#############################################
# BOLD Information
#############################################

if (! -e logs) then
     mkdir -p logs
endif
set LF = $subject_dir/$subject/logs/CBIG_preproc_create_mask.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[MASK]: logfile = $LF"
echo "Create brain mask" >> $LF
echo "[CMD]: CBIG_preproc_create_mask.csh $cmdline"   >>$LF

set boldfolder = "$subject_dir/$subject/bold"
echo "[MASK]: boldfolder = $boldfolder" |& tee -a $LF
echo "[MASK]: zpdbold = $zpdbold" |& tee -a $LF

# use the first run to create masks
cd $boldfolder

set best_run_file = "$subject_dir/$subject/qc/CBIG_preproc_bbregister_best_run.dat"
if (-e $best_run_file) then
	set mask_bold = `cat ${best_run_file}`
else
	set mask_bold = $zpdbold[1]
endif
echo "[MASK]: mask_bold = $mask_bold" |& tee -a $LF

if (! -e mask) then
	mkdir -p mask
endif

#############################################
# Create fMRI mask (wb,wm,gm,csf) according to user option
#############################################


set boldfile = $subject"_bld"$mask_bold$MASK_stem
echo "[MASK]: boldfile = $boldfile" |& tee -a $LF
set reg = $subject"_bld"$mask_bold$REG_stem".dat"
echo "[MASK]: reg = $reg" |& tee -a $LF

if ( "$erode_space" == "func" ) then
	if ( $wm == 1 && "$wm_max_erode" == "" ) then
		set wm_max_erode = 1
		echo "WARNING: WM mask will be created but -wm_max_erode is not passed in. Maximal erosion will be set to 1. (Erode in functional space.)" |& tee -a $LF
	endif
	
	if ( $csf == 1 && "$csf_max_erode" == "" ) then
		set csf_max_erode = 0
		echo "WARNING: Ventricles mask will be created but -csf_max_erode is not passed in. Maximal erosion will be set to 0. (Erode in functional space.)" |& tee -a $LF
	endif
	
else
	if ( $wm == 1 && "$wm_max_erode" == "" ) then
		set wm_max_erode = 3
		echo "WARNING: WM mask will be created but -wm_max_erode is not passed in. Maximal erosion will be set to 3. (Erode in anatomical space.)" |& tee -a $LF
	endif
	
	if ( $csf == 1 && "$csf_max_erode" == "" ) then
		set csf_max_erode = 1
		echo "WARNING: Ventricles mask will be created but -csf_max_erode is not passed in. Maximal erosion will be set to 1. (Erode in anatomical space.)" |& tee -a $LF
	endif
endif
	
# create whole brain mask	
if( $whole_brain == 1 ) then
	echo "=======================Create whole brain mask=======================" |& tee -a $LF
	if( (! -e mask/$subject.brainmask.bin.nii.gz) || ( $force == 1 ) ) then
		set cmd = "mri_vol2vol --reg $mask_bold/$reg --targ $anat_dir/$anat/mri/brainmask.mgz --mov $mask_bold/$boldfile.nii.gz --inv --o mask/$subject.brainmask.nii.gz"
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF

		set cmd = "mri_binarize --i mask/$subject.brainmask.nii.gz --o mask/$subject.brainmask.bin.nii.gz --min .0001"
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		
		echo "[MASK]: The mask file is mask/$subject.brainmask.bin.nii.gz" |& tee -a $LF
	else
		echo "[MASK]: The mask file mask/$subject.brainmask.bin.nii.gz already exists" |& tee -a $LF
	endif
endif

# create loose whole brain mask
if( $loose_whole_brain == 1 ) then
	echo "=======================Create loose whole brain mask=======================" |& tee -a $LF
	if( (! -e mask/$subject.loosebrainmask.bin.nii.gz) || ( $force == 1 ) ) then
		set cmd = "mri_vol2vol --reg $mask_bold/$reg --targ $anat_dir/$anat/mri/brainmask.mgz "
		set cmd = "$cmd --mov $mask_bold/$boldfile.nii.gz --inv --o mask/$subject.loosebrainmask.nii.gz"
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF

		set cmd = "mri_binarize --i mask/$subject.loosebrainmask.nii.gz --o mask/$subject.loosebrainmask.bin.nii.gz "
		set cmd = "$cmd --dilate 2 --min .0001"
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		
		echo "[MASK]: The mask file is mask/$subject.loosebrainmask.bin.nii.gz" |& tee -a $LF
	else
		echo "[MASK]: The mask file mask/$subject.loosebrainmask.bin.nii.gz already exists" |& tee -a $LF
	endif
endif


# create wm mask
if( $wm == 1 ) then
	echo "======================= Create wm mask =======================" |& tee -a $LF
	if( (! -e mask/$subject.func.wm.nii.gz) || ( $force == 1 )) then
		if ( "$erode_space" == "func" ) then
			set cmd = "mri_label2vol --seg $anat_dir/$anat/mri/aparc+aseg.mgz --temp $mask_bold/$boldfile.nii.gz --reg $mask_bold/$reg --o mask/$subject.func.aseg.nii" 
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
			
			foreach i (`seq 0 1 ${wm_max_erode}`)
				set cmd = "mri_binarize --i mask/$subject.func.aseg.nii --wm --erode ${i} --o mask/$subject.func.wm_erode${i}.nii.gz"
				echo $cmd |& tee -a $LF
				eval $cmd |& tee -a $LF
				
				set wm_vol = `fslstats mask/${subject}.func.wm_erode${i}.nii.gz -V`
				echo $wm_vol[2] > mask/${subject}_func_wm_erode${i}_vol.txt
			end
			
			echo "Final white matter mask is eroded by 1 times in functional space."
			set cmd = "cp mask/$subject.func.wm_erode1.nii.gz mask/$subject.func.wm.nii.gz"
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
			
		else
			foreach i (`seq 0 1 ${wm_max_erode}`)
				set cmd = "mri_binarize --i $anat_dir/$anat/mri/aparc+aseg.mgz --wm --erode ${i} --o mask/${subject}.anat.wm_erode${i}.nii.gz"
				echo $cmd |& tee -a $LF
				eval $cmd |& tee -a $LF
				
				set cmd = "mri_label2vol --seg mask/${subject}.anat.wm_erode${i}.nii.gz --temp $mask_bold/$boldfile.nii.gz --reg $mask_bold/$reg --o mask/${subject}.func.wm_erode${i}.nii.gz"
				echo $cmd |& tee -a $LF
				eval $cmd |& tee -a $LF
				
				set wm_vol = `fslstats mask/${subject}.func.wm_erode${i}.nii.gz -V`
				set wm_vol = $wm_vol[2]
				echo $wm_vol > mask/${subject}_func_wm_erode${i}_vol.txt
				set cmp_result = `echo "$wm_vol > 2700" | bc`
				if ( $cmp_result == 1 ) then
					set wm_erode = "$i"
				endif
				unset cmp_result
			end
			
			echo "Final white matter mask is eroded by $wm_erode times in anatomical space." |& tee -a $LF
			set cmd = "cp mask/$subject.func.wm_erode${wm_erode}.nii.gz mask/$subject.func.wm.nii.gz"
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
		endif
		
		echo "[MASK]: The wm mask file is mask/$subject.func.wm.nii.gz"	|& tee -a $LF
	else
		echo "[MASK]: The wm mask file mask/$subject.func.wm.nii.gz already exists" |& tee -a $LF
	endif
endif

# create csf mask
if( $csf == 1 ) then
	echo "======================= Create csf mask =======================" |& tee -a $LF
	if( (! -e mask/$subject.func.ventricles.nii.gz) || ( $force == 1 )) then
		if ( "$erode_space" == "func" ) then
			set cmd = "mri_label2vol --seg $anat_dir/$anat/mri/aparc+aseg.mgz --temp $mask_bold/$boldfile.nii.gz --reg $mask_bold/$reg --o mask/$subject.func.aseg.nii" 
			echo $cmd |& tee -a $LF
			eval $cmd
			
			set vent_erode = 0
			foreach i (`seq 0 1 ${csf_max_erode}`)  # When erode_space="func", default csf_max_erode=0. 
			                                        # Actually $subject.func.ventricles_erode1.nii.gz, ... won't be generated.
				set cmd = "mri_binarize --i mask/$subject.func.aseg.nii --ventricles --erode ${i} --o mask/$subject.func.ventricles_erode${i}.nii.gz"
				echo $cmd |& tee -a $LF
				eval $cmd
				
				set vent_vol = `fslstats mask/${subject}.func.ventricles_erode${i}.nii.gz -V`
				set vent_vol = "$vent_vol[2]"
				echo $vent_vol > mask/${subject}_func_ventricles_erode${i}_vol.txt
			end
			
			echo "Final ventricles mask is eroded by 0 times in functional space." |& tee -a $LF
			set cmd = "cp mask/$subject.func.ventricles_erode0.nii.gz mask/$subject.func.ventricles.nii.gz"
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
		
		else
			set vent_erode = 0
			foreach i (`seq 0 1 $csf_max_erode`)
				set cmd = "mri_binarize --i $anat_dir/$anat/mri/aparc+aseg.mgz --ventricles --erode ${i} --o mask/${subject}.anat.ventricles_erode${i}.nii.gz"
				echo $cmd |& tee -a $LF
				eval $cmd |& tee -a $LF
				
				set cmd = "mri_label2vol --seg mask/${subject}.anat.ventricles_erode${i}.nii.gz --temp $mask_bold/$boldfile.nii.gz --reg $mask_bold/$reg --o mask/${subject}.func.ventricles_erode${i}.nii.gz"
				echo $cmd |& tee -a $LF
				eval $cmd |& tee -a $LF
				
				set vent_vol = `fslstats mask/${subject}.func.ventricles_erode${i}.nii.gz -V`
				set vent_vol = "$vent_vol[2]"
				echo $vent_vol > mask/${subject}_func_ventricles_erode${i}_vol.txt
				set cmp_result = `echo "$vent_vol > 2700" | bc`
				if ( $cmp_result == 1 ) then
					set vent_erode = "$i"
				endif
				unset cmp_result
			end
			
			echo "Final ventricles mask is eroded by $vent_erode times in anatomical space." |& tee -a $LF
			set cmd = "cp mask/$subject.func.ventricles_erode${vent_erode}.nii.gz mask/$subject.func.ventricles.nii.gz"
			echo $cmd |& tee -a $LF
			eval $cmd |& tee -a $LF
		endif
		
		echo "[MASK]: The ventricles mask file is mask/$subject.func.ventricles.nii.gz"	|& tee -a $LF
	else
		echo "[MASK]: The ventricles mask file mask/$subject.func.ventricles.nii.gz already exists" |& tee -a $LF
	endif
endif

# create gm mask
if( $gm == 1 ) then
	echo "======================= Create gm mask=======================" |& tee -a $LF
	if( (! -e mask/$subject.func.gm.nii.gz) || ( $force == 1 )) then
		set cmd = "mri_label2vol --seg $anat_dir/$anat/mri/aparc+aseg.mgz --temp $mask_bold/$boldfile.nii.gz --reg $mask_bold/$reg --o mask/$subject.func.aseg.nii" 
		echo $cmd |& tee -a $LF
		eval $cmd

#gm labels: https://github.com/neurodebian/freesurfer/blob/cc12577ed0a75115d645a55291ae61f0d937dd7e/mri_binarize/mri_binarize.c
#match table: /apps/arch/Linux_x86_64/freesurfer/5.3.0/FreeSurferColorLUT.txt
		set cmd = "mri_binarize --i mask/$subject.func.aseg.nii --match 2 41 77 251 252 253 254 255 7 46 4 5 14 43 44 15 72 31 63 0 24 --inv --o mask/$subject.func.gm.nii.gz"
		echo $cmd |& tee -a $LF
		eval $cmd
		
		echo "[MASK]: The gm mask file is mask/$subject.func.gm.nii.gz"	|& tee -a $LF
	else
		echo "[MASK]: The gm mask file mask/$subject.func.gm.nii.gz already exists" |& tee -a $LF
	endif
endif

# If aCompCor == 1. combine wm and csf mask
if ( $aCompCor == 1 ) then
	echo "======================================= Combine wm and ventricles masks  ========================================" |& tee -a $LF
	set wm_vent_mask = "mask/${subject}.func.wm.vent.nii.gz"
	if ( ! -e ${wm_vent_mask} ) then
		set cmd = "fslmaths mask/${subject}.func.wm.nii.gz -max mask/${subject}.func.ventricles.nii.gz ${wm_vent_mask}"
		echo $cmd |& tee -a $LF
		eval $cmd |& tee -a $LF
		
		echo "The wm + ventricles mask is ${wm_vent_mask}" |& tee -a $LF
	else
		echo "The wm + ventricles mask ${wm_vent_mask} already exists" |& tee -a $LF
	endif
	
	echo "================================== aCompCor_mask done =======================================" |& tee -a $LF
endif

echo "======================= Mask done!=======================" |& tee -a $LF
echo "" |& tee -a $LF




#########################
# Output last commit of current function 
#########################
# check if git exists
which git
if (! $status) then
	echo "======================= Git: Last Commit of Current Function =======================" |& tee -a $LF
	git -C ${CBIG_CODE_DIR} log -1 -- stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_create_mask.csh >> $LF
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
			
		#anatomical name
		case "-anat_s":
			if ($#argv == 0) goto arg1err;
			set anat = $argv[1]; shift;
			breaksw
			
		#anatomical path
		case "-anat_d":
			if ($#argv == 0) goto arg1err;
			set anat_dir = $argv[1]; shift;
			breaksw
			
		#bold run number	
		case "-bld":
			if ( $#argv == 0 ) goto arg1err;
			set zpdbold = ($argv[1]); shift;
			breaksw
			
		#registration file stem
		case "-REG_stem":
			if ( $#argv == 0 ) goto arg1err;
			set REG_stem = $argv[1]; shift;
			breaksw
			
		#input file stem
		case "-MASK_stem":
			if ( $#argv == 0 ) goto arg1err;
			set MASK_stem = $argv[1]; shift;
			breaksw
		
		#create whole brain mask
		case "-whole_brain":
			set whole_brain = 1;
			breaksw	

		#create loose whole brain mask
		case "-loose_whole_brain":
			set loose_whole_brain = 1;
			breaksw
			
		#create wm mask
		case "-wm":
			set wm = 1; 
			breaksw
		
		# The maximal times to erode white matter mask
		case "-wm_max_erode":
			if ( $#argv == 0 ) goto arg1err;
			set wm_max_erode = $argv[1]; shift;
			breaksw
			
		#create csf mask
		case "-csf":
			set csf = 1; 
			breaksw
			
		# The maximal times to erode ventricles mask
		case "-csf_max_erode":
			if ( $#argv == 0 ) goto arg1err;
			set csf_max_erode = $argv[1]; shift;
			breaksw
			
		# Whether erosion is done in anatomical space or functional space
		case "-erode_space":
			if ( $#argv == 0 ) goto arg1err;
			set erode_space = $argv[1]; shift;
			breaksw
			
		#create gm mask	
		case "-gm":
			set gm = 1;
			breaksw
			
		# If -aCompCor is passed in, combine wm and ventricles masks
		case "-aCompCor":
			set aCompCor = 1;
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
if ( "$REG_stem" == "" ) then
	echo "ERROR: REG stem not specified"
	exit 1;
endif	
if ( "$MASK_stem" == "" ) then
	echo "ERROR: MASK stem not specified"
	exit 1;
endif
if ( "$erode_space" != "func" && "$erode_space" != "anat" ) then
	echo "ERROR: Wrong input for -erode_space. Please input func or anat."
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
	CBIG_preproc_create_mask.csh
	
DESCRIPTION:
	This function use T1-T2* registration file to create fMRI masks (wb,wm,gm,csf).
	
	The masks are created based on FreeSurfer segmentation. 

	The whole brain mask is created by mri_vol2vol using anatomical volume 
	<anat_dir>/<anat_src>/mri/brainmask.mgz as target. It is binarized by a 
	threshold of 0.0001. The loose whole brain mask is created by dilating 
	the whole brain mask 2 voxels.

	The white matter mask is created from aseg 2 and 41. The ventricles mask 
	is created from aseg ventricles + choroid (not include 4). The grey matter 
	is created by excluding aseg 2 41 77 251 252 253 254 255 7 46 4 5 14 43 44 
	15 72 31 63 0 24.
	
	If the WM mask and ventricles mask are eroded in functional space, the apac+aseg 
	will be first transformed into functional space by mri_label2vol, then the 
	corresponding segments are chosen and eroded by mri_binarize. 
	If the WM mask and ventricles mask are eroded in anatomical space, the 
	corresponding segments in anatomical apac+aseg will be chosen and eroded by 
	mri_binarize, then it will be transformed into functional space by mri_label2vol.
	
	If the WM mask and ventricles mask are eroded in functional space, the final 
	mask for WM is eroded once, and the final mask for ventricles is not eroded.
	If the WM mask and ventricles mask are eroded in anatomical space, the final 
	mask for WM is the smallest (depends on the maximal erosions the user passed in) 
	but has at least 100 voxels, and the final mask for ventricles is also the one 
	which is smallest (depends on the maximal erosions the user passed in) but has 
	at least 100 voxels.

REQUIRED ARGUMENTS:
	-s  subject_id                 : name of the subject
	-d  subject_dir                : absolute path to <subject_id>, i.e. all preprocessed data of this 
	                                 subject is assumed to be stored in <subject_dir>/<subject_id>.
	-anat_s  anat_src              : name of anatomical folder for the subject, i.e. recon-all results
	-anat_d  anat_dir              : absolute path to <anat_src>, i.e. all recon-all results of this 
	                                 subject is assumed to be stored in <anat_dir>/<anat_src>.
	-bld  bold_runs                : all bold run numbers of this subject. Each number should be three 
	                                 digits. If this subject has multiple runs, please use space as 
	                                 delimiter between two runs (e.g. -bld '002 003'). NOTE: quote sign 
	                                 is necessary.
	-REG_stem  REG_stem            : stem of T1-T2* registration file. E.g. if the registration file is 
	                                 Sub0001_Ses1_bld002_rest_skip4_stc_mc_reg.dat, then <REG_stem> = 
	                                 _rest_skip4_stc_mc_reg. This registration file is assumed to be 
	                                 stored in <subject_dir>/<subject_id>/bold/<run_number>
	-MASK_stem  MASK_stem          : stem of input file used to create the mask. E.g. if the template to 
	                                 create the masks is Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz, then 
	                                 <MASK_stem> = _rest_skip4_stc_mc. Notice that the choice of <MASK_stem> 
	                                 and <REG_stem> should be reasonable and consistent. For instance, the 
	                                 user cannot apply the registration file before motion correction on 
	                                 the template volume after motion correction.

OPTIONAL ARGUMENTS:
	-whole_brain                   : if this option is used, create whole brain mask
	-loose_whole_brain             : if this option is used, create loose whole brain mask by dilating 
	                                 whole brain mask 2 voxels
	-wm                            : if this option is used, create white matter mask
	-wm_max_erode  wm_max_erode    : the maximal voxels to be eroded for WM mask. If erode in functional 
	                                 space, default is 1; if erode in anatomical space, default is 3.
	                                 All WM masks eroded from 0 to <wm_max_erode> will be created.
	-csf                           : if this option is used, create ventricle mask
	-csf_max_erode  csf_max_erode  : the maximal voxels to be eroded for ventricles mask. If erode in
	                                 functional space, default is 0; if erode in anatomical space, default 
	                                 is 1. All ventricles masks eroded from 0 to <csf_max_erode> will be
	                                 created.
	-gm                            : if this option is used, create grey matter mask
	-aCompCor                      : if this option is used, the WM mask and ventricles mask will be merged.
	-erode_space  erode_space      : The space that the erosion of WM and ventricles masks will be done. 
	                                 Choose from "func" or "anat". Default is anat.
	-force                         : update results, if exist then overwrite
	-help                          : help
	-version                       : version

OUTPUTS:
	A folder <subject_dir>/<subject_id>/bold/mask containing all the masks that the users specified to create.

EXAMPLE:
	$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_mask.csh -s Sub0001_Ses1 
	-d ~/storage/fMRI_preprocess -anat_s Sub0001_Ses1_FS 	-anat_d ~/storage/sMRI_preprocess -bld "002 003" 
	-REG_stem _rest_skip4_stc_mc_reg -MASK_stem _rest_skip4_stc_mc -whole_brain -wm -csf -gm

