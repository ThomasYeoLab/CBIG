#!/bin/tcsh -f
#
# This function runs CBIG_preproc_deoblique.sh before recon-all for input T1 (and T2 if used). Input files will be 
# changed to plumb and in RPI orientation (3dinfo) before running recon-all. Please use this function unless you are 
# sure your input data is already plumb and in RPI orientation. 
# To check your input files, use 3dinfo <input_file>, and check Data Axes Tilt and Data Axes Approximate Orientation.
# Note that mri_info will use opposite notation for orientation. RPI in 3dinfo is equivalent to LAS in mri_info.
#
# Options are the same as recon-all. You can simply replace recon-all with CBIG_preproc_recon-all.csh
# If you are not familiar with recon-all, please refer to https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all
#
# Example: CBIG_preproc_recon-all.csh -all -i T1.nii.gz -s Sub001 -sd <sub_dir>
#
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Print help and version
set VERSION = '$Id: CBIG_preproc_recon-all.csh v 1.0 2023/02/09'

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

set opt = "";

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

# Deoblique and reorient input T1 file
echo "Checking input T1 file..."
if ( $input =~ "*.nii.gz" ) then
    set output = ${input:r:r}
    set ext = "nii.gz"
else if ( $input =~ "*.nii" ) then
    set output = ${input:r}
    set ext = "nii"
endif
set output = "${output}_deoblique.${ext}"
set cmd = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/utilities/CBIG_preproc_deoblique.sh"
set cmd = "${cmd} -i ${input} -o ${output}"
echo ${cmd}
eval ${cmd}
echo "Deoblique finished for input T1."
set opt = "${opt} -i ${output}"


# Deoblique and reorient input T2 file
if ( $?input_T2 ) then
    echo "Checking input T2 file..."
    if ( $input_T2 =~ "*.nii.gz" ) then
        set output_T2 = ${input_T2:r:r}
        set ext_T2 = "nii.gz"
    else if ( $input_T2 =~ "*.nii" ) then
        set output_T2 = ${input_T2:r}
        set ext_T2 = "nii"
    endif
    set output_T2 = "${output_T2}_deoblique.${ext_T2}"
    set cmd = "$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/utilities/CBIG_preproc_deoblique.sh"
    set cmd = "${cmd} -i ${input_T2} -o ${output_T2}"
    echo ${cmd}
    eval ${cmd}
    echo "Deoblique finished for input T2."
    set opt = "${opt} -T2 ${output_T2}"
endif

# recon-all
echo "Start recon-all..."
set cmd = "recon-all ${opt}"
echo ${cmd}
eval ${cmd}

exit 0;

##################################
# Parse Arguments
##################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
    set flag = $argv[1]; shift;
    switch($flag)
        case "-i":
            if ( $#argv == 0 ) goto arg1err;
            set input = $argv[1]; shift;
            breaksw

        case "-T2":
            if ( $#argv == 0 ) goto arg1err;
            set input_T2 = $argv[1]; shift;
            breaksw    
        
        default:
            set opt = "$opt $flag";
            breaksw
    endsw
end

goto parse_args_return;

############################################
# Check Parameters
############################################
check_params:
if ( "$input" == "" ) then
    echo "ERROR: T1 is not specified."
    exit 1;
endif

goto check_params_return;


#######################################################
# ERROR message
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
# Help
#######################################################
BEGINHELP

NAME:
    CBIG_preproc_recon-all.csh

DESCRIPTION:
    Call CBIG_preproc_deoblique.sh to deoblique and reorient T1 input (and T2 if passed in) before running recon-all.

REQUIRED ARGUMENTS:
    -s          subject   : subject id
    -i          input     : input T1 nifti file

OPTIONAL ARGUMENTS:
    -T2         input_T2  : input T2 nifti file

Options are the same as recon-all. You can simply replace recon-all with CBIG_preproc_recon-all.csh
Please refer to recon-all for other optional arguments and outputs. 
If you are not familiar with recon-all, please refer to https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all

Example:
    CBIG_preproc_recon-all.csh -all -i T1.nii.gz -s Sub001 -sd <sub_dir>
