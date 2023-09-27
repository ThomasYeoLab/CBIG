#!/bin/csh -f

# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles.csh v 1.0 2016/06/18 $'

set sub_dir = ""
set sub = ""
set surf_data = ""
set roi = fsaverage3
set target = fsaverage5
set threshold = 0.1
set scrub_flag = 0;

set PrintHelp = 0;
set n = `echo $argv | grep -e -help | wc -l`
if( $#argv == 0 || $n != 0 ) then
    echo $VERSION
    # print help
    cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
    exit 0;
endif
set n = `echo $argv | grep -e -version | wc -l`
if( $n != 0 ) then
    echo $VERSION
    exit 0;
endif

set code_dir = `pwd`;

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set root_dir = `python -c "import os; print(os.path.realpath('$0'))"`
set root_dir = `dirname $root_dir`

set MATLAB = `which $CBIG_MATLAB_DIR/bin/matlab`
if($status) then
    echo "ERROR: could not find matlab"
    exit 1
endif

if( $scrub_flag == 0 ) then
    set output_file1 = "${output_dir}/lh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.nii.gz"
    set output_file2 = "${output_dir}/rh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.nii.gz"
else
    set output_file1 = "${output_dir}/lh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile_scrub.nii.gz"
    set output_file2 = "${output_dir}/rh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile_scrub.nii.gz"
endif

echo "===>> Compute correlation profile for ${sub_dir}/${sub}, writing outputs in $output_dir"

############################
# create input files
############################
set lh_input_file = "${output_dir}/lh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.input"
if( -e ${lh_input_file} ) then
    rm  $lh_input_file
endif
foreach surf ($surf_data)
    set lh_surf = `echo $surf | sed "s/rh/lh/g"`
    echo $lh_surf >> $lh_input_file
end

set rh_input_file = "${output_dir}/rh.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.input"
if( -e ${rh_input_file} ) then
    rm  $rh_input_file
endif
foreach surf ($surf_data)
    set rh_surf = `echo $surf | sed "s/lh/rh/g"`
    echo $rh_surf >> $rh_input_file
end

if( $scrub_flag == 1 ) then
    set outlier_input_file = "${output_dir}/outlier.${sub}.roi${roi}.thres${threshold}.surf2surf_profile.input"
    if( -e ${outlier_input_file} ) then
        rm  $outlier_input_file
    endif
    foreach outlier ($outlier_files)
        echo $outlier >> $outlier_input_file
    end
    echo "outlier_input_file = $outlier_input_file"
endif


###########################
# compute correlation profile
###########################

if( -e ${output_file1} && -e ${output_file2} ) then
    echo "Outputs already exist. Skipping......"
else
    set cmd = ( $MATLAB -nodesktop -nodisplay -nosplash -r '"')
    set cmd = ($cmd 'addpath(genpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"utilities"'", "'"matlab"'"')));')
    set cmd = ($cmd 'addpath(fullfile('"'"$CBIG_CODE_DIR"'", "'"external_packages"'",)
    set cmd = ($cmd "'"SD"'", "'"SDv1.5.1-svn593"'", "'"BasicTools"'"'));')
    if( $scrub_flag == 0 ) then
        set cmd = ($cmd CBIG_ComputeCorrelationProfile "'"${roi}"'" "'"${target}"'" "'"${output_file1}"'")
        set cmd = ($cmd "'"${output_file2}"'" "'"${threshold}"'" "'"${lh_input_file}"'")
        set cmd = ($cmd "'"${rh_input_file}"'"'; exit;''"')
    endif
    
    if( $scrub_flag == 1 ) then
        set cmd = ($cmd CBIG_ComputeCorrelationProfile "'"${roi}"'" "'"${target}"'" "'"${output_file1}"'")
        set cmd = ($cmd "'"${output_file2}"'" "'"${threshold}"'" "'"${lh_input_file}"'")
        set cmd = ($cmd "'"${rh_input_file}"'" "'"${outlier_input_file}"'"'; exit;''"')
    endif
    eval $cmd
endif

if( -e ${output_file1} && -e ${output_file2} ) then
    echo "===>> Matlab succeeds to compute correlation profiles."
else
    echo "===>> Matlab fails to compute correlation profiles."
endif

echo ""
exit 0;



############################
# parse arguments
############################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
    set flag = $argv[1]; shift;
    
    switch($flag)
        # subject directory
        case "-sd":
            if( $#argv == 0 ) goto arg1err;
            set sub_dir = $argv[1]; shift;
            breaksw
        
        # subject
        case "-s":
            if( $#argv == 0 ) goto arg1err;
            set sub = $argv[1]; shift;
            breaksw

        # surface data
        case "-surf_data":
            if( $#argv == 0 ) goto argerr;
            set surf_data = "$argv[1]"; shift;
            breaksw
            
        # outlier files
        case "-outlier_files":
            if( $#argv == 0 ) goto arg1err;
            set outlier_files = "$argv[1]"; shift;
            set scrub_flag = 1
            breaksw
            
        # target resolution
        case "-target":
            if( $#argv == 0 ) goto arg1err;
            set target = $argv[1]; shift;
            breaksw
            
        # ROI resolution
        case "-roi":
            if( $#argv == 0 ) goto arg1err;
            set roi = $argv[1]; shift;
            breaksw

        # output directory
        case "-output_dir":
            if( $#argv == 0 ) goto arg1err;
            set output_dir = $argv[1]; shift;
            breaksw

        default:
            echo "ERROR: Flag $flag unrecognized."
            echo $cmdline
            exit 1;
            breaksw

    endsw 
end

goto parse_args_return;


############################
# check parameters
############################
check_params:

if( $#sub_dir == 0 ) then
    echo "ERROR: subject directory not specified."
    exit 1;
endif

if( $#sub == 0 ) then
    echo "ERROR: subject id not specified."
    exit 1;
endif

if( $#surf_data == 0 ) then
    echo "ERROR: surface data not specified."
    exit 1;
endif

if( $#output_dir == 0 ) then
    echo "ERROR: output directory not specified."
    exit 1;
endif

goto check_params_return;


############################
# Error message
############################
arg1err:
    echo "ERROR: flag $flag requires one argument"
    exit 1;

argerr:
    echo "ERROR: flag $flag requires at lease one arguments"
    exit 1;


    

exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

NAME:
    CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles.csh

DESCRIPTION:
    This function computes the functional connectivity profiles on surface for the given subject. 
    The resolution of the functional connectivity is given by the inputs -taget and -roi. Default 
    is fsaverage5 x fsaverage3.
    
    The motion outliers file is optional. If motion outliers file is passed in, the high motion 
    frames will be ignored when computing functional connectivity. High-motion frames are 
    indicated by 0 in motion outliers files, while low-motion frames are 1.
  
REQUIRED ARGUMENTS:
    -sd             sub_dir       : fMRI subjects directory. This directory contains all the folders
                                    named by the subject IDs.
    -s              subject       : subject id
    -surf_data      surf_data     : surface fMRI data filenames of all runs of <subject> (only left  
                                    hemisphere, or only right hemisphere). Please use space as the 
                                    delimiter between different runs. For example, 
                                    -surf_data "<sub_dir>/<subject>/surf/lh.Sub0001_Ses1_bld002*.nii.gz \
                                    <sub_dir>/<subject>/surf/lh.Sub0001_Ses1_bld003*.nii.gz"
                                    NOTE: quote sign is necessary.
                                    This function will transform all "lh" to "rh" when computing left 
                                    hemisphere correlation, or "rh" to "lh" when computing right 
                                    hemisphere correlation.
    -output_dir     output_dir    : directory to output FC profiles file (full path)

OPTIONAL ARGUMENTS:
    -outlier_files  outlier_files : motion outliers files of all runs. Please use space as the delimiter 
                                    between different runs. For example, 
                                    -outlier_files "<sub_dir>/<subject>/qc/lh.Sub0001_Ses1_bld002*.txt \
                                    <sub_dir>/<subject>/qc/lh.Sub0001_Ses1_bld003*.txt"
                                    NOTE: quote sign is necessary.
    -target         target        : the resolution of clustering (default is fsaverage5)
    -roi            roi           : the resolution of ROIs (default is fsaverage3)

OUTPUTS:
    "<output_dir>/lh.<subject>.roifsaverage3.thres0.1.surf2surf_profile.input" 
    is a text file where each line is the surface fMRI data (lh) of one run of <subject>.
    
    "<output_dir>/rh.<subject>.roifsaverage3.thres0.1.surf2surf_profile.input" 
    is a text file where each line is the surface fMRI data (rh) of one run of <subject>.
    
    "<output_dir>/outlier.<subject>.roifsaverage3.thres0.1.surf2surf_profile.input"
    is a text file where each line is the motion outlier filename of one run of <subject>.
    
    "<output_dir>/lh.<subject>.roifsaverage3.thres0.1.surf2surf_profile_scrub.nii.gz"
    is the surface to surface correlation file (lh) of <subject>. This correlation profile
    is computed by choosing a seed in the mesh of <target> from left hemisphere and correlated
    with all ROIs in the mesh of <roi> from both left and right hemisphere.
    
    "<output_dir>/rh.<subject>.roifsaverage3.thres0.1.surf2surf_profile_scrub.nii.gz"
    is the surface to surface correlation file (rh) of <subject>. This correlation profile
    is computed by choosing a seed in the mesh of <target> from right hemisphere and correlated
    with all ROIs in the mesh of <roi> from both left and right hemisphere.

EXAMPLE:
    CBIG_compute_fcMRI_surf2surf_profiles.csh -sd ~/storage/fMRI_data -s Sub0001_Ses1 
    -surf_data "~/storage/fMRI_data/Sub0001_Ses1/surf/ \
    lh.Sub0001_Ses1_bld002_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz 
    ~/storage/fMRI_data/Sub0001_Ses1/surf/ \
    lh.Sub0001_Ses1_bld003_rest_skip4_stc_mc_resid_cen_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz" 
    -outlier_files "~/storage/fMRI_data/Sub0001_Ses1/qc/Sub0001_Ses1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt 
    ~/storage/fMRI_data/Sub0001_Ses1/qc/Sub0001_Ses1_bld003_FDRMS0.2_DVARS50_motion_outliers.txt"
    -target fsaverage5 -roi fsaverage3

