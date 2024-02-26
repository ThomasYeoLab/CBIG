#!/bin/csh -f

###############################################
# # convert BIDS for fMRI preprocessing
# #############################################
# #############################################
# # In this script, we: 
# # 1. find relevant .nii(.gz) and .json files from BIDS subject folder
# # 2. check consistency between image header and json files
# # 3. create .fmrinii file which can be used for CBIG_preprocessing pipeline
# #
# # Example:  
# # $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_convert_input_BIDS.csh
# # -d ~/storage/test_dataset -s A001 -stem "task-rest" -d "input_conversion_pipeline"
# #############################################
# # Written by Yan Hongwei and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_convert_input_BIDS.sch, v 1.0 2023/09/25'

set n = `echo $argv | grep -e -help |wc -l`

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

set dataset_dir = ""     # full path to BIDS dataset
set subject_id = ""      # BIDS subject ID (without prefix sub-)
set session = ""         # keyword to filter session (with prefix ses-)
set stem = ""            # keywords to filter images. (e.g. "run-1 rest")
set output_dir = "CBIG_preproc_convert_input_BIDS" 
                         # output of this function will be saved here
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

#####################################################
# Create log files
#####################################################


# if output_dir is relative path, create output_dir according to BIDS convention
if (`echo $output_dir | cut -c1-1` != "/") then
    set output_dir = "$dataset_dir/derivatives/$output_dir/"
endif

set log_dir = "$output_dir/sub-$subject_id/logs"
if ( "$session" != "" ) then
    set log_dir = $output_dir"/sub-"$subject_id"_"$session"/logs"
endif
if (! -e $log_dir ) then
    mkdir -p $log_dir
endif

set LF = $log_dir/CBIG_preproc_convert_input_BIDS.log
if ( -e $LF ) then
    rm $LF
endif
touch $LF
echo "[convert_input_BIDS]: logfile = $LF"
echo "convert_input_BIDS" >> $LF
echo "[CMD]: CBIG_preproc_convert_input_BIDS.csh $cmdline" >> $LF

set fmrinii_dir = "$output_dir/fmrinii"
if ( ! -e $fmrinii_dir) then
    mkdir -p $fmrinii_dir
endif

#####################################################
# find .nii(.gz) files
#####################################################

echo ">>> Searching for .nii(.gz) files in subject folder>>>" >> $LF
cd $dataset_dir/sub-$subject_id

# only files in func/ directory is relevant
set image_list = `find "$PWD" \( -name '*_bold.nii.gz' -o -name '*_bold.nii' \) | grep '\/func\/' | sort`

# filter with session, if it is provided
if ("$session" != "") then
    set image_list = `echo $image_list | fmt -1 | grep $session`
endif

# filter with stem, if it is provided
if ( "$stem" != "") then
    set keyword_list = `echo $stem:q | sed 's/ / /g'`
    foreach keyword ($keyword_list)
        set image_list = `echo $image_list | fmt -1 | grep $keyword`
    end
endif

set number_of_image = $#image_list
if ( $number_of_image == 0 ) then
    echo "[ERROR]: No relevant image found. Exit." | & tee -a $LF
    exit 1
endif

echo "Relevant images files: " >> $LF
@ i = 1
while ( $i <= $number_of_image )
    echo "$image_list[$i]" >> $LF
    @ i++
end

####################################################
# check consistency between json and image files
####################################################

echo ">>>Check consistency between JSON files and image header>>>" >> $LF
@ i = 1
while ( $i <= $number_of_image )
    # search for matching json files, starting from func/ level and moves up to dataset/ level
    # at func/ level, json file should have the same filename as image file, except for extension
    # at other levels, json filename should end with "_bold.json"
    
    # func/ level 
    set found_json = 0
    if ( "$image_list[$i]" =~ *.nii.gz ) then
        set json_file = `echo $image_list[$i] | sed 's/.\{6\}$//'`"json"    #".nii.gz"
    else
        set json_file = `echo $image_list[$i] | sed 's/.\{3\}$//'`"json"    #".nii"
    endif
    if ( -e $json_file ) then
        set found_json = 1
    else
        # other levels
        set is_dataset_dir = 0
        set curr_dir = `dirname $image_list[$i]`
        cd $curr_dir
        while ( $is_dataset_dir == 0 && $found_json == 0 )
            cd ..
            set json_file_list = `find "$PWD" -maxdepth 1 -name '*_bold.json'`
            if ( $#json_file_list > 1 ) then
                echo "[WARNING]: More than one relevant json file found for $image_list[$i]" | & tee -a $LF
            endif
            if ( $#json_file_list > 0 ) then
                set json_file = $json_file_list[1]
                set found_json = 1
            endif
            set curr_dir = $PWD
            if ( `realpath $curr_dir` == `realpath $dataset_dir`) then
                set is_dataset_dir = 1
            endif
        end
    endif

    if ( $found_json == 1 ) then
        echo "Identified $json_file as relevant JSON file for $image_list[$i], checking consistency.." >> $LF
        # check TR consistency
        set TR_json = `cat $json_file | tr -d '\n' | awk -F '[:,}]' '/"RepetitionTime"/{print $2}' | tr -d '"'`
        set TR_json_convert = `echo " $TR_json * 1000 / 1 " | bc`
        set TR_header = `mri_info --tr $image_list[$i]`
        if ( ( $TR_json != $TR_header ) && ( $TR_json_convert != $TR_header ) ) then
            set msg = "[WARNING]: TR from $json_file is not consistent with TR from $image_list[$i] header."
            set msg = "$msg TR from json: $TR_json, TR from image header: $TR_header."
            set msg = "$msg  Will follow TR from image header."
            echo "$msg" | & tee -a $LF
        else
            echo "JSON file and image header are consistent." >> $LF
        endif
        # future note: add more checks when needed, json files will be ignored in the pipeline
    else
        echo "No relevant JSON file found for $image_list[$i], skip JSON consistency check" >> $LF
    endif
    @ i++
end

#####################################################
# check if there is only one session
#####################################################

echo ">>> Checking number of sessions in detected files>>>" >> $LF

# note: /ses-*/ folder is optional for BIDS dataset
# i.e. folder structure can be either
# <dataset>/sub-<subject_id>/func/ or
# <dataset>/sub-<subject_id>/ses-*/func/

# search for "ses-*" keyword in paths
set ses_list = ()
@ i = 1
while ( $i <= $number_of_image )
    set token_list = `echo $image_list[$i] | sed 's/\// /g'`
    set filtered_token_list = ()
    foreach token ($token_list)
        if ( `echo $token | grep 'ses-' | grep -v '_'` != "") then
            set filtered_token_list = ($filtered_token_list $token)
        endif
    end
    if ( $#filtered_token_list > 1 ) then
        echo "[WARNING]: Multiple session folders detected in file path: $image_list[$i]" >> $LF
    endif
    if ( $#filtered_token_list >= 1 ) then
        set ses_list = ($ses_list $filtered_token_list[1])
    endif
    @ i++
end
set ses_list = `echo $ses_list | fmt -1 | sort | uniq`

if ( $#ses_list > 1 ) then
    echo "[ERROR]: Multiple sessions detected in relevant image list." | & tee -a $LF
    if ("$session" == "") then
        echo "[ERROR]: Add the argument '-ses' to select one session." | & tee -a $LF
    else
        echo "[ERROR]: Modify the argument '-ses' to select one session." | & tee -a $LF
    endif
    echo "[ERROR]: See log file for detail."
    exit 1
endif

#####################################################
# create .fmrinii file
#####################################################

echo ">>>Generate .fmrinii file>>>" >> $LF
if ( "$session" != "" ) then
    set session = "_$session"
endif
set fmrinii_file = $fmrinii_dir"/sub-"$subject_id$session".fmrinii"
if ( -e $fmrinii_file ) then
    rm $fmrinii_file
endif
touch $fmrinii_file
echo ".fmrinii filename: $fmrinii_file" >> $LF

# check if there are multiple runs for the subject
# Note: the entity 'run-' might not exist in filename

# make a copy of image_list and remove "echo-*" from the filename
set copied_image_list = ()
@ i = 1
while ( $i <= $number_of_image )
    set token_list = `echo $image_list[$i] | sed 's/_/ /g'`
    set copied_image = ""
    foreach token ($token_list)
        if ( "`echo $token | grep -c '^echo-'`" == "0" ) then
            set copied_image = "$copied_image$token"
        endif
    end
    set copied_image_list = ($copied_image_list $copied_image)
    @ i++
end

# find unique names in copied_image_list and record its indexes
set unique_name_list = ($copied_image_list[1])
set index_list = (1) #index of each unique name in original list
@ index_count = 1  #total number of unique names
@ i = 2
while ( $i <= $number_of_image )
    set is_found = 0
    @ j = 1
    while ( $j <= $#unique_name_list )
        if ( $copied_image_list[$i] == $unique_name_list[$j] ) then
            set index_list = ($index_list $j)
            set is_found = 1
        endif
        @ j++
    end
    if ( $is_found == 0 ) then
        set unique_name_list = ($unique_name_list $copied_image_list[$i])
        @ index_count++
        set index_list = ($index_list $index_count)
    endif
    @ i++
end

# add file paths to .fmrinii file
if ( $index_count == 1 ) then
    # no session folder or only one session is relevant
    echo "assuming all images are from a single run." >> $LF
    echo "001 $image_list" >> $fmrinii_file
else
    # multiple runs
    echo "multiple runs detected, adding multiple rows to .fmrinii" >> $LF
    @ i = 1
    while ( $i <= $index_count )
        set run_image_list = ()
        @ j = 1
        while ( $j <= $number_of_image )
            if ($index_list[$j] == $i) then
                set run_image_list = ($run_image_list $image_list[$j])
            endif
            @ j++
        end
        echo "`printf "%03.0f" $i` $run_image_list" >> $fmrinii_file
        @ i++
    end
endif

echo "fmrinii file generated!" | & tee -a $LF

exit 0;

#####################################################
# Parse Arguments
#####################################################

parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
    set flag = $argv[1]; shift;

    switch($flag)
        # dataset dir
        case "-d":
            if ( $#argv == 0 ) goto arg1err;
            set dataset_dir = $argv[1]; shift;
            breaksw

        # subject name
        case "-s":
            if ( $#argv == 0 ) goto arg1err;
            set subject_id = $argv[1]; shift;
            breaksw

        # output directory
        case "-output_dir":
            if ( $#argv == 0 ) goto arg1err;
            set output_dir = $argv[1]; shift;
            breaksw

        # session of subject
        case "-ses":
            if ( $#argv == 0 ) goto arg1err;
            set session = $argv[1]; shift;
            breaksw

        # stem to filter filenames in subject folder
        case "-stem":
            if ( $#argv == 0 ) goto arg1err;
            set stem = $argv[1]; shift;
            breaksw

        default:
            echo "ERROR: Flag $flag unrecognized."
            echo $cmdline
            exit 1
            breaksw
    endsw
end

goto parse_args_return;

#####################################################
# Check Parameters
#####################################################

check_params:
if ( "$dataset_dir" == "" ) then
    echo "ERROR: dataset directory is not specified"
    exit 1
endif
if ( "$subject_id" == "" ) then
    echo "ERROR: BIDS subject not specified"
    exit 1
endif
# make sure subject_id does not start with "sub-"
if ( "$subject_id" =~ sub-* ) then
    set correct_subject_id = `echo $subject_id | cut -c 5-`
    set msg = "WARNING: BIDS subject_id should not start with 'sub-'."
    set msg = "$msg Assuming subject_id is $correct_subject_id and subject folder is $subject_id"
    echo $msg
    set subject_id = $correct_subject_id
endif
# make sure session starts with "ses-"
if ( "$session" != "" && "`echo $session | grep -c '^ses-'`" == "0" ) then
    set session = "ses-$session"
endif
goto check_params_return;

#####################################################
# ERROR message
#####################################################

arg1err:
    echo "ERROR: flag $flag requires one argument"
    exit 1

#####################################################
# Help
#####################################################
BEGINHELP

NAME:
        CBIG_preproc_convert_input_BIDS.csh

DESCRIPTION:
        Convertor for BIDS formatted dataset.
        This function does:
        1) find relevant .nii(.gz) and .json files from BIDS subject folder
        2) check consistency between image header and json files
        3) create .fmrinii file which can be used for CBIG_preprocessing pipeline. The .fmrinii file includes 
           all _bold.nii(.gz) files under <dataset>/sub-<subject_id>/*/func/ directories. If <stem> is specified,
           the files will be filtered such that only files matching all keywords are included. Each line in 
           the .fmrinii file corresponds to one run. If there are multiple echos in the run, the files for the 
           same run will appear in the same line, separated by one white space.
        IMPORTANT: If there are multiple sessions in the subject folder, then -ses must be provided. If -ses
        is not provided, this function assumes there is only one session in the subject folder.

REQUIRED ARGUMENTS:
        -d <dataset_dir>         : absolute path to <dataset_dir>. This directory is assumed to be BIDS
                                   compliant
        -s <subject_id>          : subject ID (without 'sub-' prefix). fMRI data for the subject is assumed to 
                                   be at <dataset_dir>/sub-<subject_id>/*/func/

OPTIONAL ARGUMENTS:
        -ses <session>           : <session> is used to filter the session. The preprocessing pipeline processes
                                   one session each time. If there is only one session folder in the subject folder,
                                   or if there is no session folder in the subject folder, this argument is not 
                                   required.
                                   The argument can be with or without prefix 'ses-'
        -stem <stem>             : <stem> is used to filter nii files in sub_<subject_id> directory.
                                   If not specified, all .nii.gz) files will be used. If more than one keywords
                                   are required, please use space as delimiter between keywords 
                                   (e.g. -stem "run-1 rest") NOTE: quote sign is necessary. 
        -output_dir <output_dir> : default value is 'CBIG_preproc_convert_input_BIDS'.
                                   If <output_dir> is a relative path, then the output of this script will be 
                                   saved in <dataset_dir>/derivatives/<output_dir>. 
                                   If <output_dir> is an absolute path, then the output of this script will be
                                   saved in <output_dir>.
                                   In <output_dir>, fmrinii/ directory will be created for the .fmrinii file,
                                   sub-<subject_id>/logs/ directory will be created for log file.
                                   It is recommended to use default value or relative path for BIDS compliance.

OUTPUTS:
        This function will produce the following files:
        1) <output_dir>/fmrinii/sub-<subject_id>(_ses-<session>).fmrinii This file can be used as input for 
           CBIG_preproc_fMRI_preprocess.csh
        2) <output_dir>/sub-<subject_id>(_ses-<session>)/logs/CBIG_preproc_convert_input_BIDS.log This file contains
           log for this function.

EXAMPLE:
        $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_convert_input_BIDS.csh 
-d ~/storage/test_dataset -s A001 -ses 1 -stem "task-rest" -d "input_conversion_pipeline"
