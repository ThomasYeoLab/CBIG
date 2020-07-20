#!/bin/csh -f

#############################################
# Spatial distortion correction with fieldmap
#############################################
#############################################
# In this script, we:
# 0) Assume the functional image has already been motion corrected and transformation matrices from
#    motion correction are available
# 1) Process the fieldmap based on the given image
#    The fieldmap can have one of the following two forms:
#    a. a phase difference image and a magnitude image
#    b. two opposite phase encoding direction fieldmap images (e.g: AP and PA)
# 2) Unwarp the functional image using FUGUE:
#    This step applies the warping field generated from the processed fieldmaps to 
#    the functional image
# The unwarping method is taken from Midnight Scan Club (MSC) preprocessing pipeline, epi_unwarp_MSC,
# which is written by Tyler Blazey 
# (https://github.com/MidnightScanClub/MSCcodebase/blob/master/Processing/epi_unwarp_MSC).
#############################################
# Author: Shaoshi Zhang
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_spatial_distortion_correction.csh, v 1.0 2018/11/1 $'

set n = `echo $argv | grep -e -help | wc -l`

# if there is no arguments or there is -help option 
if( $#argv == 0 || $n != 0 ) then
        echo $VERSION
        # print help    
        cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
        exit 0;
endif

set subject = ""
set sub_dir = ""
set zpdbold = ""
set BOLD_stem = ""

set fpm = ""                #fieldmap preocessing method (mag+phasediff/opposite PED)
set mag = ""                #magnitude image absoulte path (for mag+phasediff)
set phase = ""              #phase-differene image absolute path (for mag+phasediff)
set delta_te = ""           #difference of te between two fieldmap magnitude image (for mag+phasediff)

set j_plus = ""             #j+ direction image absolute path (for oppo PED)
set j_minus = ""            #j- direction image absolute path (for oppo PED)
set j_minus_trt = ""        #total readout time of j- image (for oppo PED, unit is in seconds)
set j_plus_trt = ""         #toral readout time of j+ image (for oppo PED, unit is in seconds)

set ees = ""                #effective echo spacing of fMRI image (for both, unit is in microseconds)
set TE = ""                 #echo time of fMRI image (for both, unit is in microseconds)

set ref = 0                 #reference frame, default is the first frame
set prelude = 1             #unwrap the phase-difference image
#do not change unit when processing the fieldmap if set to 1(will be set to 1 if fieldmap processing method is oppo PED)
set no_unit_conversion = 0	
set topup_config = ""       #config file for topup (use default if not specified)
set sig_threshold = 0.9     #signal loss threshold. Default is 0.9
set fmri_bet = 0.2          #BET fractional intensity threshold for fMRI image. Default is 0.2
set mag_bet = 0.3           #BET fractional intensity threshold for magnitude image. Default is 0.3



		
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

cd $sub_dir/$subject

#############################################
# BOLD Information
#############################################

set qc = $sub_dir/$subject"/qc"

if (! -e $qc) then
     mkdir -p $qc
endif
if (! -e logs) then
     mkdir -p logs
endif
set LF = $sub_dir/$subject/logs/CBIG_preproc_spatial_distortion_correction.log
if( -e $LF ) then
	rm $LF
endif
touch $LF
echo "[SDC]: logfile = $LF"
echo "Spatial Distortion Correction" >> $LF
echo "[CMD]: CBIG_preproc_spatial_distortion_correction.csh $cmdline"   >>$LF

set boldfolder = "$sub_dir/$subject/bold"
set sdc = "$boldfolder"/sdc
if ( ! -e $sdc ) then
     mkdir $sdc
endif
echo "[SDC]: boldfolder = $boldfolder" |& tee -a $LF
echo "[SDC]: zpdbold = $zpdbold" |& tee -a $LF

#############################################
# Process fieldmap
#############################################
echo "=======================Processing input fieldmap images=======================" |& tee -a $LF

cd $sdc
#############################################
# When fieldmaps are in magnitude and phasediff
#############################################
if ( $fpm == "mag+phasediff" ) then
     if ( $mag == "" ) then
          echo "ERROR: fieldmap magnitude image not specified. Exit." |& tee -a $LF
          exit 1
     endif
     if ( $phase == "" ) then
          echo "ERROR: fieldmap phase-difference image not specified. Exit." |& tee -a $LF
          exit 1
     endif
     if ( $delta_te == "" ) then
          echo "ERROR: delta TE not specified. Exit." |& tee -a $LF
          exit 1
     endif
     imcp $mag $sdc/magnitude
     imcp $phase $sdc/phase_difference
     set mag = magnitude
     set phase = phase_difference

     # Create brain mask using magnitude image
     set cmd = (bet $mag "$mag"_brain -f $mag_bet -m -R)
     echo $cmd |& tee -a $LF
     eval $cmd 
     set mag = "$mag"_brain
     set mag_brain_mask = "$mag"_mask

	# Unit converstion to rad/s (and unwrapping) 
	if ( $no_unit_conversion == 0 ) then
		set range = `fslstats ${phase} -R`
		set min = `echo $range[1] | cut -d '.' -f1`
		set max = `echo $range[2] | cut -d '.' -f1`
		set range = `echo $min $max`
		if ( $range[1] >= -4096 && $range[2] <= 4094 ) then
			set cmd = ( fslmaths $phase -div 4096 -mul 3.14159265 "$phase"_rad -odt float )
		else if ( $range[1] == 0 && $range[2] == 4094 ) then
			set cmd = ( fslmaths $phase -sub 2047 -div 2047 -mul 3.14159265 "$phase"_rad )
			set cmd = ( $cmd -odt float )
		else if ( $range[1] == 0 && $range[2] == 4095 ) then
			set cmd = ( fslmaths $phase -sub 2047.5 -div 2047.5 -mul 3.14159265 "$phase"_rad )
			set cmd = ( $cmd -odt float )
		else
			echo "Phase image has values that are not within expected range. Exit." |& tee -a $LF
			exit 1
		endif
		echo $cmd |& tee -a $LF
		eval $cmd
		set phase = "$phase"_rad

		#Unwrap phase image if prelude is set
		if ( $prelude == 1 ) then
			set cmd = ( prelude -p $phase -a $mag -o "$phase"_unwrap -m "$mag"_mask )
			echo $cmd |& tee -a $LF
			eval $cmd
			set phase = "$phase"_unwrap
		endif

		#divide by delta TE 
		set cmd = ( fslmaths $phase -div "$delta_te" "$phase"_sec )
		echo $cmd |& tee -a $LF
		eval $cmd
		set phase = "$phase"_sec
	endif
	echo "Fieldmap magnitude image			$mag" |& tee -a $LF
	echo "Fieldmap phase image				$phase" |& tee -a $LF
	echo "Fieldmap magnitude brain mask		$mag_brain_mask" |& tee -a $LF

#############################################
# When two fieldmaps have opposite phase encoding directions
#############################################
else if ( $fpm == "oppo_PED" ) then
	if ( $j_plus == "" ) then
		echo "ERROR: j+ image not specified. Exit."|& tee -a $LF
		exit 1
	endif
	if ( $j_minus == "" ) then
		echo "ERROR: j- image not specified. Exit."|& tee -a $LF
		exit 1
	endif
	if ( $j_minus_trt == "" || $j_plus_trt == "" ) then
		echo "ERROR: Total readout time not specified. Exit. "|& tee -a $LF
		exit 1
	endif
	#if user does not specify topup config file, use default config file
	if ( $topup_config == "" ) then
		echo "[WARNING] FSL TOPUP configuration file not specified. Use default topup config file."|& tee -a $LF
		set topup_config = "$CBIG_FSLDIR"/src/topup/flirtsch/b02b0.cnf
	endif

	#create datain.txt
	set datain = "$sdc"/datain.txt
	if ( -e $datain ) then
	    rm $datain 
	endif
	touch $datain
	echo "0 -1 0 $j_minus_trt" >> $datain	#AP/j-
	echo "0 1 0 $j_plus_trt" >> $datain	#PA/j+
	echo "[SDC] datain.txt successfully created." |& tee -a $LF

	#merge two opposite PED images (j- followed by j+)
	set cmd = ( fslmerge -t AP_PA.nii.gz "$j_minus" "$j_plus" )
	echo $cmd |& tee -a $LF
	eval $cmd
	#run TOPUP
	set cmd = ( topup --imain=AP_PA.nii.gz --datain="$datain" --config="$topup_config" --fout=fmap_phase_hz.nii.gz  )
	set cmd = ($cmd --iout=fmap_unwarped.nii.gz )
	echo $cmd |& tee -a $LF
	eval $cmd
	#convert the unit to rad/s
	set cmd = ( fslmaths fmap_phase_hz -mul 6.28 fmap_phase_rad_sec )
	echo $cmd |& tee -a $LF
	eval $cmd
	set phase = fmap_phase_rad_sec
	#average across magnitude images and brain extraction
	set cmd = ( fslmaths fmap_unwarped -Tmean fmap_mag )
	echo $cmd |& tee -a $LF
	eval $cmd
	set cmd = ( bet fmap_mag fmap_mag_brain -f $mag_bet -m -R )
	echo $cmd |& tee -a $LF
	eval $cmd
	set mag = fmap_mag_brain
	set mag_brain_mask = fmap_mag_brain_mask

	echo "Fieldmap magnitude image 		$mag" |& tee -a $LF
	echo "Fieldmap phase image	 		$phase" |& tee -a $LF
	echo "Fieldmap magnitude brain mask	$mag_brain_mask" |& tee -a $LF

# If the fieldmap proceesing method is neither 'mag+phasediff' nor 'oppo_PED', exit with error
else 
	echo "ERROR: Unknown fieldmap processing method. (fieldmap processing method should be either 'mag+phasediff' \
	or 'oppo_PED'). Exit."|& tee -a $LF
	exit 1
endif

echo "=======================Fieldmap processing done=======================" |& tee -a $LF
echo "" |& tee -a $LF



##########################################
# Unwarp BOLD images
########################################## 
foreach curr_bold ($zpdbold) 
     cd $boldfolder/$curr_bold
     set boldfile = "$subject"_bld"$curr_bold""$BOLD_stem"
#copy over processed magnitude, phase fieldmap image, and brain mask
     set cmd = ( imcp $sdc/$mag_brain_mask fmap_magnitude_brain_mask )
     eval $cmd
     set cmd = ( imcp $sdc/$mag fmap_magnitude_brain )
     eval $cmd
     set cmd = ( imcp $sdc/$phase fmap_phase )
     eval $cmd
end

foreach curr_bold ($zpdbold) 
cd $boldfolder/$curr_bold
set boldfile = "$subject"_bld"$curr_bold""$BOLD_stem"

echo "=======================Create Brain Mask=======================" |& tee -a $LF

set mag_brain_mask = fmap_magnitude_brain_mask
set mag = fmap_magnitude_brain
set phase = fmap_phase

#Extract the reference frame
echo "Reference frame for $boldfile.nii.gz --> $ref" |& tee -a $LF
set cmd = ( fslroi $boldfile "$boldfile"_"$ref" $ref 1 )
echo $cmd |& tee -a $LF
eval $cmd

#fMRI brain extraction
set cmd = ( bet "$boldfile"_"$ref" "$boldfile"_"$ref"_brain -m -f $fmri_bet -R )
echo $cmd |& tee -a $LF
eval $cmd

#Create inverted field map brain mask
set cmd = ( fslmaths $phase -abs -bin -mas "$mag_brain_mask" -mul -1 -add 1 -bin "$phase"_inv_brain_mask )
echo $cmd |& tee -a $LF
eval $cmd
		 
#Cluster inverted brain mask
set cmd = ( cluster -i "$phase"_inv_brain_mask -t 0.5 --no_table -o "$phase"_inv_brain_mask_clust )
echo $cmd |& tee -a $LF
eval $cmd 
		
#Save intensity of largest cluster 
set max = `fslstats "$phase"_inv_brain_mask_clust -R | awk '{print $2}'`

#Threshhold the image by max, then invert again. Get a new, tighter brain mask.
set cmd = ( fslmaths "$phase"_inv_brain_mask_clust -thr $max -bin -mul -1 -add 1 -bin )
set cmd = ( $cmd -mas "$mag_brain_mask" "$mag_brain_mask" )
echo $cmd |& tee -a $LF
eval $cmd 
		 
#Use the new brain mask on the phase image
set cmd = ( fslmaths "$phase" -mas "$mag_brain_mask" "$phase"_masked )
echo $cmd |& tee -a $LF
eval $cmd
set phase = "$phase"_masked 

#Get a 50% brain mask
set fifty = `fslstats "$mag" -P 98 | awk '{print ( $1 / 2 ) }'`
set cmd = ( fslmaths "$mag" -thr $fifty -bin "$mag"_fifty_mask )
echo $cmd |& tee -a $LF
eval $cmd 

#Erode the original brain mask
set cmd = ( fslmaths "$mag_brain_mask" -ero "$mag_brain_mask"_eroded )
echo $cmd |& tee -a $LF
eval $cmd 

#Add eroded and fifty masks
set cmd = ( fslmaths "$mag_brain_mask"_eroded -add "$mag"_fifty_mask -thr 0.5 -bin "$mag_brain_mask" )
echo $cmd |& tee -a $LF
eval $cmd 
		 
#Mask the phase image again
set cmd = ( fslmaths "$phase" -mas "$mag_brain_mask" "$phase" )
echo $cmd |& tee -a $LF
eval $cmd 

#Erode brain mask again
set cmd = ( fslmaths "$mag_brain_mask" -ero "$mag_brain_mask"_eroded )
echo $cmd |& tee -a $LF
eval $cmd 

echo "=======================Brain Mask Created=======================" |& tee -a $LF
echo "" |& tee -a $LF
echo "=======================Apply filetring on phase image=======================" |& tee -a $LF

#Create filter
set filter = "$phase"_filter_despike
set cmd = ( fugue --loadfmap="$phase" --mask="$mag_brain_mask" --despike --despikethreshold=2.1 --savefmap="$filter" )
echo $cmd |& tee -a $LF
eval $cmd 

#Apply the filter to brain edges
set cmd = ( fslmaths "$filter" -sub "$phase" -mas "$mag_brain_mask"_eroded -add "$phase" "$phase"_filtered )
echo $cmd |& tee -a $LF
eval $cmd 
set phase = "$phase"_filtered

#Shift median to 0
set median = `fslstats "$phase" -k "$mag_brain_mask" -P 50`
set cmd = ( fslmaths "$phase" -sub $median "$phase"_norm )
echo $cmd |& tee -a $LF
eval $cmd
set phase = "$phase"_norm
echo "=======================Filtering done=======================" |& tee -a $LF
echo "" |& tee -a $LF

echo "=======================Calculating warping field=======================" |& tee -a $LF
set xfm = "$boldfile"_"$ref"_brain_to_fmap_mag_brain_signal_lossed_distorted.mat
set inv_xfm = fmap_mag_brain_signal_lossed_distorted_to_"$boldfile"_"$ref"_brain.mat
	
#Estimate signal loss from phase image. Range goes from 0 (no signal) to 1 (full signal)
set cmd = ( sigloss -i "$phase" --te="$TE" -m "$mag_brain_mask" -s "$phase"_signal_loss )
echo $cmd |& tee -a $LF 
eval $cmd 

#Multiply fieldmap magnitude image by signal loss image. Will result in a magnitude with areas of signal loss.
set cmd = ( fslmaths "$phase"_signal_loss -mul "$mag" fmap_mag_brain_signal_lossed -odt float )
echo $cmd |& tee -a $LF
eval $cmd 

#Run fugue on the signal lossed magnitude image. Will distort it according to the fieldmap phase image. 
set cmd = ( fugue -i fmap_mag_brain_signal_lossed --loadfmap="$phase" )
set cmd = ( $cmd --mask="$mag_brain_mask" -w fmap_mag_brain_signal_lossed_distorted )
set cmd = ( $cmd --nokspace --unwarpdir='y-' --dwell="$ees" )
echo $cmd |& tee -a $LF
eval $cmd 

#Do the same thing for the signal_loss phase image.
set cmd = ( fugue -i "$phase"_signal_loss --loadfmap="$phase" --dwell="$ees" )
set cmd = ( $cmd -w "$phase"_signal_loss_distorted --nokspace --unwarpdir='y-' )
set cmd = ( $cmd --mask="$mag_brain_mask" )
echo $cmd |& tee -a $LF
eval $cmd 

#Threshhold the distorted signal loss brain according to user chosen signal loss threshhold
set cmd = ( fslmaths "$phase"_signal_loss_distorted -thr "$sig_threshold" )
set cmd = ( $cmd "$phase"_signal_loss_distorted )
echo $cmd |& tee -a $LF
eval $cmd 

#Register the distorted magnitude to the distorted functional image.
#Generate func --> fmap_mag 
#Use the threshholded signal loss brain as weighting (areas with 0 will be ignored)
set cmd = ( flirt -ref fmap_mag_brain_signal_lossed_distorted -in "$boldfile"_"$ref"_brain )
set cmd = ( $cmd -omat $xfm -refweight "$phase"_signal_loss_distorted -dof 6 )
set cmd = ( $cmd -out "$boldfile"_"$ref"_brain_to_fmap_mag_brain_signal_lossed_distorted )
echo $cmd |& tee -a $LF
eval $cmd 

#Invert transformation to get fmap_mag --> func
set cmd = ( convert_xfm -omat $inv_xfm -inverse $xfm )
echo $cmd |& tee -a $LF
eval $cmd 

#Apply transformation to original magnitude image for QA purposes
set cmd = ( flirt -in "$mag" -ref "$boldfile"_"$ref"_brain -applyxfm -init )
set cmd = ( $cmd $inv_xfm -out fmap_mag_brain_to_"$boldfile"_"$ref"_brain )
echo $cmd |& tee -a $LF
eval $cmd 

#Apply transformation (fmap_mag --> func) to the field phase image, the result will be 
#phase image registered to functional image
set cmd = ( flirt -in "$phase" -ref "$boldfile"_"$ref"_brain -applyxfm )
set cmd = ( $cmd -init $inv_xfm -out "$phase"_to_"$boldfile"_"$ref"_brain )
echo $cmd |& tee -a $LF
eval $cmd 
set phase = "$phase"_to_"$boldfile"_"$ref"_brain

#Apply transform to mask
set cmd = ( flirt -in "$mag_brain_mask" -ref "$boldfile"_"$ref"_brain -applyxfm )
set cmd = ( $cmd -init $inv_xfm -out fmap_mag_brain_mask_to_"$boldfile"_"$ref"_brain )
echo $cmd |& tee -a $LF
eval $cmd  

#Rebinarize mask
set cmd = ( fslmaths fmap_mag_brain_mask_to_"$boldfile"_"$ref"_brain -thr 0.5 -bin )
set cmd = ( $cmd fmap_mag_brain_mask_to_"$boldfile"_"$ref"_brain )
echo $cmd |& tee -a $LF
eval $cmd |& tee 
			  
#dilate brain mask slightly in order to prevent erosion
set cmd = ( fslmaths fmap_mag_brain_mask_to_"$boldfile"_"$ref"_brain -dilM )
set cmd = ( $cmd -dilM -ero fmap_mag_brain_mask_to_"$boldfile"_"$ref"_brain )
echo $cmd |& tee -a $LF
eval $cmd  

#Run fugue using the registerted phase image on the distorted func. 
#Will unwarp the func reference frame and save the shift map that does this.
set cmd = ( fugue --loadfmap="$phase" --dwell=$ees -u "$boldfile"_"$ref"_unwarped )
set cmd = ( $cmd  -i "$boldfile"_"$ref" --saveshift="$boldfile"_"$ref"_unwarp_shift )
set cmd = ( $cmd --mask=fmap_mag_brain_mask_to_"$boldfile"_"$ref"_brain --unwarpdir='y-' )
echo $cmd |& tee -a $LF
eval $cmd 

#Convert the shiftwarp to an absolute warp
set cmd = ( convertwarp -s "$boldfile"_"$ref"_unwarp_shift -r "$boldfile"_"$ref" )
set cmd = ( $cmd --shiftdir='y-' -o "$boldfile"_"$ref"_unwarp_shift_warp ) 
echo $cmd |& tee -a $LF
eval $cmd 

echo "=======================Calculate warping field finished=======================" |& tee -a $LF
echo "" |& tee -a $LF
echo "=======================Unwarp BOLD image=======================" |& tee -a $LF
#check existence of relevant input images for applywarp
if ( ! -e "$boldfile"_"$ref"_unwarp_shift_warp.nii.gz ) then
    echo "[ERROR] Warping field missing, please check if FSL version is 5.0.10 or above." |& tee -a $LF
    exit 1
endif

if ( ! -e "$boldfile"_mc.cat ) then
    echo "[ERROR] motion correction transformation matrices missing!" |& tee -a $LF
    exit 1
endif
#apply motion correction transformation matrices and unwarping field at the same time to
#reduce the number of interpolation

set cmd = ( applywarp -i "$boldfile" -o "$boldfile"_mc_sdc -r "$boldfile" --abs )
set cmd = ( $cmd -w "$boldfile"_"$ref"_unwarp_shift_warp --premat="$boldfile"_mc.cat --interp=spline )
echo $cmd |& tee -a $LF
eval $cmd

echo "=======================Unwarping done!=======================" |& tee -a $LF
echo "" |& tee -a $LF

#tidy up the space
if ( ! -e unwarping ) then
     mkdir warping
endif
mv fmap* warping/
mv "$boldfile"_"$ref"* warping/

end

exit 0

##########################################
# Parse Arguments 
##########################################      

parse_args:
set cmdline = "$argv";

while( $#argv != 0 )
        set flag = $argv[1]; shift;

        switch($flag)
		#subject id
		case "-s":
			if ( $#argv == 0 ) goto argerr;
			set subject = $argv[1]; shift;
			breaksw

		#subject folder
		case "-d":
			if ( $#argv == 0 ) goto argerr;
			set sub_dir = $argv[1]; shift;
			breaksw

		#bold run number
		case "-bld":
			if ( $#argv == 0 ) goto argerr;
			set zpdbold = ($argv[1]); shift;
			breaksw
		
		#input file stem
		case "-BOLD_stem":
			if ( $#argv == 0 ) goto argerr;
			set BOLD_stem = $argv[1]; shift;
			breaksw
		
		#fieldmap processing method
		case "-fpm":
			if ( $#argv == 0 ) goto argerr;
			set fpm = $argv[1]; shift;
			breaksw

		#magnitude image absolute path
		case "-m":
			set mag = $argv[1]; shift;
			breaksw

		#phase-difference image absolute path
		case "-p":
			set phase = $argv[1]; shift;
			breaksw

		#delta TE
		case "-delta":
			set delta_te = `echo "$argv[1]/1000" | bc -l`; shift;
			breaksw
		
		#j+ image absolute path
		case "-j_plus":
			set j_plus = $argv[1]; shift;
			breaksw

		#j- image absolute path
		case "-j_minus":
			set j_minus = $argv[1]; shift;
			breaksw

		#total readout time j-
		case "-j_minus_trt":
			set j_minus_trt = $argv[1]; shift;
			breaksw

		#total readout time j-
		case "-j_plus_trt":
			set j_plus_trt = $argv[1]; shift;
			breaksw
		
		#effective echo spacing
		case "-ees":
			if ( $#argv == 0 ) goto argerr;
			set ees = `echo "$argv[1]/1000" | bc -l`; shift;
			breaksw

		#echo time
		case "-te":
			if ( $#argv == 0 ) goto argerr;
			set TE = `echo "$argv[1]/1000" | bc -l`; shift;
			breaksw

		#reference frame
		case "-ref":
			set ref = $argv[1]; shift;
			breaksw

		#use prelude to unwrap phase image
		case "-no_prelude":
			set prelude = 0; shift;
			breaksw
		
		#do not convert fieldmap unit (i.e. already in rad/s)
		case "-no_unit_conversion":
			set no_unit_conversion = 1; shift;
			breaksw

		#topup configuraion file absolute path
		case "-topup_config":
			set topup_config = $argv[1]; shift;
			breaksw

		#signal loss threshold
		case "-sig":
			set sig_threshold = $argv[1]; shift;
			breaksw

		#BET threshold for fMRI image
		case "-fmri_bet":
			set fmri_bet = $argv[1]; shift;
			breaksw

		#BET threshold for magnitude image
		case "-mag_bet":
			set mag_bet = $argv[1]; shift;
			breaksw
		#cleanup intermediate files
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
if ( "$fpm" == "" ) then
	echo "ERROR: fieldmap processing method not specified"
	exit 1;
endif
if ( "$ees" == "" ) then
	echo "ERROR: effective echo spacing of fMRI image not specified"	
	exit 1;
endif
if ( "$TE" == "" ) then
	echo "ERROR: TE of fMRI image not specified"
	exit 1;
endif

goto check_params_return;

##########################################
# ERROR message
##########################################      

argerr:
  echo "ERROR: flag $flag requires one argument"
  exit 1

#####################################
# Help
#####################################
BEGINHELP

NAME:
	CBIG_preproc_spatial_distortion_correction
DESCRIPTION:
	 This function:
  	 0) assumes the fMRI image has already been motion corrected and transformation matrices available
	 1) processes the filedmap based on the given image
		1.1) If the input images are in the form of magnitude and phase-difference, this function will
		     create a mask based on the magnitude image, and convert the unit of the phase-differene image 
		     to rad/s. Unwrapping the phase image is optional.
		1.2) If the input images are along opposite phase encoding direction, this function will run FSL TOPUP
		     to unwarp the images and genreate an esimated field in HZ, the unit will be converted to rad/s. 
		     Currently, only phase encoding directions along 'j-' and 'j' are supported. Ensure FSL version is at least 
		     5.0.10, otherwise the output fieldmap image may have wrong orientation. For now, this script is only
			 tested on FSL 5.0.10.
			 [NOTE] Check the orientation of fieldmap to determine whether the phase encoding direction is 'j-' or 'j'.
			 For example, if the orientation of fieldmap follows RAS convention, this means voxel position increases from
			 posterior to anterior, therefore a phase encoding direction 'PA' is equivalent to 'j' and vice versa.
	 2) Calculates warping field of the fMRI image using FUGUE
		First, this function creates a brain mask based on the fieldmap magnitude image and applies it on the phase image.
		Second, the phase image is further filetered by despiking, demeaned, and forward warped using FUGUE.
		Third, this functions register phase image to fucntional image by the following step
		Register the functional image to magnitude image and obtain the transformation
		Inverse the above transformation and get the registration from magnitude image to functional image
		Then apply the above transformation to phase image, such that phase image is registered to functional image
		Fourth, run FUGUE using the registered phase image to obtain the warping field
	 3) Unwarps the fMRI image
	    fMRI image is unwarpped using FSL applywarp, this step requires the motion correction 
	    transformation matrices beforehand. The matrices should be stored in a single file with name 
	    <sub>_bld<bold_run>_<bold_stem>_mc.cat under <sub_dir>/<sub/>bold/<bold_run>.
	    Warping field and motion correction transformation matrices will be applied together to 
	    reduce the number of interpolation.
	    
	 The method of calculating the warping field is taken from Midnight Scan Club (MSC) preprocessing pipeline,
	 epi_unwarp_MSC, which is written by Tyler Blazey 
	 (https://github.com/MidnightScanClub/MSCcodebase/blob/master/Processing/epi_unwarp_MSC).

	 See spatial_distortion_correction_readme for more details. 


REQUIRED ARGUMENTS:
	-s  <subject_id>           			: subject's ID
	-d  <subject_dir>         			: absolute path to <subject_id>. More specifically, processed data of 
	                             	  	  this subject will be in folder <subject_dir>/<subject_id>
	-bld  <bold_runs>          			: bold run numbers, each number must have three digits. If this 
	                            	  	  subject has multiple runs, please use space as delimiter between 
	                                  	  two run numbers (e.g. -bld "002 003"). NOTE: quote sign is necessary.
	-BOLD_stem  <BOLD_stem>   			: specify the stem of input file. (e.g. if input file is 
	                           		  	  Sub0001_bld002_rest_skip4_stc.nii.gz then the stem of this file is 
	                           			  <BOLD_stem> = _rest_skip4_stc). This input file is assumed to be 
	                           		  	  stored in <subject_dir>/<subject_id>/bold/<run_number>.
	-fpm <fieldmap_processing_method>	: fieldmap processing method, can be mag+phasediff or oppo_PED
	-ees <effective echo spacing>		: effective echo spacing of the functional image (ms)
	-te <echo time>						: echo time of the functional image (ms)

OPTIONAL ARGUMENTS:
	-m <magnitude image>				: absolute path of the magnitude image
	-p <phase image>					: absolute path of the phase image
	-delta <delta_te>					: differene of echo time of two magnitude images (ms) 
	-j_minus <j- image>					: absolute path of the j- image
	-j_plus <j+ image>					: absolute paht of the j image
	-j_minus_trt <total readout time>	: total readout time of j- image (s)
	-j_plus_trt <total readout time>	: total readout time of j image (s)
	-ref <ref>							: reference frame, default is the first frame (which starts with 0)
	-no_unit_conversion					: do not perform unit conversion to rad/s
	-no_prelude							: do not unwrap the phase image
	-topup_config <config>				: absolute path of topup configuration file, use default if not specified
	-sig <signal loss threshold>		: signal loss threshold, default is 0.1
	-fmri_bet <fmri_bet>				: BET threshold for fMRI image, default is 0.2
	-mag_bet <mag_bet>					: BET threshold for fieldmap magnitude image, default is 0.3
	-help								: help
	-version							: version
	

OUTPUTS:
	(1) A NIFTI volume after spatial distortion correction
	    <sub_dir>/<subject>/bold/<bold_run>/<subject>_bld<bold_run>_<BOLD_stem>_mc_sdc.nii.gz
	(2) Unwarping folder containins:
		fmap_magnitude_brain		:brain extracted magnitude image, obtained from fieldmap and used for masking purpose
		fmap_magnitude_brain_mask	:mask obtained from fmap_magnitude_brain
		fmap_phase					:phase image obtained from fieldmap, unit is in rad/s
		fmap_phase_masked			:phase image that is masked by fmap_magnitude_brain_mask
		fmap_phase_masked_filtered_norm_signal_loss_distorted:
						phase image that is masked, despiked, demeaned, signal loss estimated, and forward warped.
						This image is further used as reference weight to register magnitude image to functional image
		<sub>_bld<bold_run>_<BOLD_stem>_ref<ref>:
						a reference frame from functional image
		<sub>_bld<bold_run>_<BOLD_stem>_ref<ref>_unwarped:
						an unwarped reference frame
		<sub>_bld<bold_run>_<BOLD_stem>_ref<ref>_unwarp_shift:
						a voxel shfit map obtained based on the unwarped frame
		<sub>_bld<bold_run>_<BOLD_stem>_ref<ref>_unwarp_shift_warp:
						the warping field used to unwarp all frames of functional image
		fmap_mag_brain_to_<sub>_bld<bld_rin>_<BOLD_stem>_ref<ref>.nii.gz:
						registration from magnitude image to functional image
		fmap_phase_masked_filetered_norm_to_<sub>_bld<bld_rin>_<BOLD_stem>_ref<ref>.nii.gz:
						registration from phase image to functional image
		
EXAMPLE:
	./CBIG_preproc_spatial_distortion_correction -s sub-NDARAA536PTU -d ~/storage/fMRI_preprocess -bld '001' \
	-BOLD_stem _rest -fpm "mag+phasediff" -m /data/HBN/rawData_release1_4/SI/sub-NDARAA536PTU/fmap/ \
	sub-NDARAA536PTU_magnitude1.nii.gz -p /data/HBN/rawData_release_1_4/SI/sub-NDARAA536PTU/fmap/ \
	sub-NDARAA536PTU_phasediff.nii.gz -delta 4.76 -ees 0.55 -te 40 

	./CBIG_preproc_spatial_distortion_correction.csh -s sub-NDARWN691CG7 -d ~/storage/fMIR_preprocess -bld '001 002'\
	 -BOLD_stem _rest -fpm oppo_PED -j_minus /data/HBN/rawData_release1_4/RU/sub-NDARWN691CG7/fmap/ \
	sub-NDARWN691CG7_dir-AP_acq-fMRI_epi.nii.gz -j_plus /data/HBN/rawData_release1_4/RU/sub-NDARWN691CG7/\
	fmap/sub-NDARWN691CG7_dir-PA_acq-fMRI_epi.nii.gz -j_minus_trt 0.04565 -j_plus_trt 0.04565 -ees .580013000 -te 30.00


