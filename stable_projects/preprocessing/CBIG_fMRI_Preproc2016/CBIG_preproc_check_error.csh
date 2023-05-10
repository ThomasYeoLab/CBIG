#!/bin/csh -f

# example:
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_check_error.csh 
#    -sd ~/storage/fMRI_preprocess -sub_ls ~/storage/fMRI_preprocess/scripts/subject_list.txt 
#    -out_list ~/storage/fMRI_preprocess/scripts/unsuccessful_subject_list.txt 
#    -out_msg ~/storage/fMRI_preprocess/scripts/unsuccessful_msg.txt
# 
# This function checks error message within all log files for a single subject.
#
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$Id: CBIG_preproc_check_error.csh v 1.0  $'

set sub_dir = ""
set sub_list = ""
set output_file = ""

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

goto parse_args;
parse_args_return:

goto check_params;
check_params_return:

set out_list_dir = `dirname $out_list`
mkdir -p $out_list_dir
if( -e ${out_list} ) then
    rm $out_list
endif
touch $out_list

set out_msg_dir = `dirname $out_msg`
mkdir -p $out_msg_dir
if( -e ${out_msg} ) then
    rm $out_msg
endif
touch $out_msg

set subjects = `cat ${sub_list}`

set s_err = ""
foreach s ($subjects)
    echo "$s :"
    if ( ! -d $sub_dir/$s ) then
        echo "$s   does not have its output folder. This subject hasn't start preprocess." |& tee -a $out_msg
        echo $s |& tee -a $out_list

    else
        set err_flag = 0
        set succ_flag = 0

        # log file directory
        set log_dir = ${sub_dir}/${s}/logs

        # all .log file inside this directory
        set logfiles = `ls ${log_dir}/*.log`
        #echo "all log files: $logfiles"

        # check each log file if there exists "ERROR" or "Maximum number of clients reached"
        set log_err = ""
        foreach log ($logfiles)
            set logbase = `basename $log`

            set find_str = `egrep -c 'ERROR|Error|Maximum number of clients reached|Segmentation violation detected' ${log}`

            if ( "$find_str" != "0" ) then
                if ( "$logbase" != "git.log" ) then
                    set err_flag = 1
                    set log_err = "${log_err}  ${logbase}"
                endif
            endif
        end

        # check main log file if there is successful message
        set main_log = "${log_dir}/CBIG_preproc_fMRI_preprocess.log"
        set find_succ = `egrep -c 'Preprocessing Completed!' ${main_log}`
        if ( "$find_succ" != "0" ) then
            set succ_flag = 1
        endif

        if ( $err_flag == 1 ) then
            if ( $succ_flag == 0 ) then
                echo "$s    $log_err  have error messages." |& tee -a $out_msg
                echo "$s" |& tee -a $out_list
            else
                echo "$s    $log_err  have error messages, but the final successful message exists." |& tee -a $out_msg
                echo "$s" |& tee -a $out_list
            endif
        else
            if ( $succ_flag == 0 ) then
                echo "$s    No error message, but terminated by unexpected reason." |& tee -a $out_msg
                echo "$s" |& tee -a $out_list
            endif
        endif
    endif
end


exit 0


#####################################
# parse arguments
#####################################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
    set flag = $argv[1]; shift;

    switch($flag)
        # subjects directory
        case "-sd":
            if( $#argv == 0 ) goto arg1err
            set sub_dir = $argv[1]; shift;
            breaksw

        # subject id
        case "-sub_ls":
            if( $#argv == 0 ) goto arg1err
            set sub_list = $argv[1]; shift;
            breaksw

        # output subject list
        case "-out_list":
            if( $#argv == 0 ) goto arg1err;
            set out_list = $argv[1]; shift;
            breaksw

        # output message file
        case "-out_msg":
            if( $#argv == 0 ) goto arg1err;
            set out_msg = $argv[1]; shift;
            breaksw

        default:
            echo "ERROR: flag $flag unrecognized."
            echo $cmdline
            exit 1
            breaksw
    endsw
end

goto parse_args_return;


########################################
# check arguments
########################################
check_params:

if( $#sub_dir == 0 ) then
    echo "ERROR: subjects directory not specified."
    exit 1;
endif

if( $#sub_list == 0 ) then
    echo "ERROR: subject if not specified."
    exit 1
endif

if( $#output_file == 0 ) then
    echo "ERROR: output file name not specified."
    exit 1
endif

goto check_params_return;


########################################
# Error message
########################################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1;
  

#######################################
# usage exit
#######################################
BEGINHELP

NAME:
    CBIG_preproc_check_error.csh

DESCRIPTION:
    Check error message in log files generated through CBIG_preproc_fMRI_preprocess preprocessing pipeline. 

    This function takes in the subjects directory and subjects list and output a text file contains 
    the information of unsuccessful subjects. 

    If one or more log files of a subject contains error messages ('ERROR' or 'Error' or 'Maximum 
    number of clients reached'), the name of this subject and the its log files with error messages 
    will be printed into the output text file. 

    Each subject with log file names occupies one line in the output text file.

REQUIRED ARGUMENTS:
    -sd        sub_dir     : absolute path to fMRI subjects' folders. E.g. the preprocessed 
                             data of subject Sub0001_Ses1 are stored in <sub_dir>/Sub0001_Ses1. 
                             The folder structure follows the one given by CBIG_preproc_fMRI_preprocess.csh
    -sub_ls    sub_list    : subjects list (full path), each line is one subject ID.
    -out_list  out_list    : output subject list (full path)
    -out_msg   out_msg     : output message file name (full path)

OUTPUTS:
    The output text files given by -out_list and -out_msg options.

EXAMPLE:
    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_check_error.csh 
    -sd ~/storage/fMRI_preprocess -sub_ls ~/storage/fMRI_preprocess/scripts/subject_list.txt 
    -out_list ~/storage/fMRI_preprocess/scripts/unsuccessful_subject_list.txt 
    -out_msg ~/storage/fMRI_preprocess/scripts/unsuccessful_msg.txt

