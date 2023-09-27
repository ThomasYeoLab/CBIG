#!/bin/csh -f

# Given the subjects directory, the subjects list, a relative data folder (e.g. bold, vol, surf), and the data stem, this function create a list (or two lists for surface data) where each line is one subject with different runs.
# This function assumes the fMRI data are preprocessed by CBIG_preproc_fMRI_preprocess.csh.
# The user needs to specify the output directory and output list stem.
# Date: 2016/07/17
# Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set VERSION = '$ID: CBIG_fMRI_create_data_list.csh, v 1.0 2016/06/18'

set sub_dir = "";
set sub_list = "";
set folder = ""
set data_stem = ""
set out_dir = ""
set out_stem = ""
set preproc_opt = "new"

set PrintHelp = 0;
if( $#argv == 0 ) goto usage_exit;
set n = `echo $argv | grep -e -help | wc -l`
if( $n != 0 ) then
    set PrintHelp = 1;
    goto usage_exit;
endif
set n = `echo $argv | grep -e -version | wc -l`
if( $n != 0 ) then
    echo $VERSION
    exit 0;
endif


goto parse_args;
parse_args_return:

goto check_params;
check_params_return:


#################################
# output configuration
#################################
mkdir -p ${out_dir}
if( $folder == "surf" ) then
    set lh_output = "${out_dir}/lh.${out_stem}.list"
    set rh_output = "${out_dir}/rh.${out_stem}.list"
    if( -e ${lh_output} ) then
        echo "WARNING: ${lh_output} already exists. Overwriting the old one."
        rm ${lh_output}
    endif
    if( -e ${rh_output} ) then
        echo "WARNING: ${rh_output} already exists. Overwriting the old one."
        rm ${rh_output}
    endif

else
    set output = "${out_dir}/${out_stem}.list"
    if( -e ${output} ) then
        echo "WARNING: ${output} already exists. Overwriting the old one."
        rm ${output}
    endif
endif

pushd $sub_dir > /dev/null

#################################
# loop through all subjects, loop through runs of each subject
#################################
set subjects = `cat $sub_list`;
foreach sub ($subjects)
    pushd $sub > /dev/null

    if( $preproc_opt == "new" ) then
        if( ! -e ./logs/${sub}.bold ) then
            echo "ERROR: bold list ${sub_dir}/${sub}/logs/${sub}.bold does not exist!"
            exit 1;
        endif
        set bold = `cat ./logs/${sub}.bold`
    endif
    if( $preproc_opt == "old" ) then
        pushd $sub_dir/$sub/scripts > /dev/null
        eval "`grep "fcbold" *.params`"

        set zpdbold = ""
        @ k = 1
        while ($k <= ${#fcbold})
            set zpdbold = ($zpdbold `echo $fcbold[$k] | awk '{printf ("%03d",$1)}'`)
            @ k++
        end

        set bold = "$zpdbold"
        popd > /dev/null
    endif

    if( $folder != "surf" && $folder != "bold") then
        pushd ${folder} > /dev/null
        set data_names = ""
        foreach bold_num ($bold)
            set curr_data_name = "${sub}_bld${bold_num}${data_stem}"
            set curr_data_name = "${sub_dir}/${sub}/${folder}/${curr_data_name}"
            if( $#curr_data_name != 1 ) then
                echo "ERROR: No data file or multiple data files with ${data_stem} for bold run number ${bold_num}. Please check file existance or the bold list."
                exit 1;
            endif

            set data_names = "${data_names} ${curr_data_name}"

        end
        echo ${data_names} >> $output
        popd > /dev/null
    endif

    if( $folder == "bold" ) then
        pushd ${folder} > /dev/null
        set data_names = ""
        foreach bold_num ($bold)
            set curr_data_name = "${sub}_bld${bold_num}${data_stem}"
            set curr_data_name = "${sub_dir}/${sub}/${folder}/${bold_num}/${curr_data_name}"
            if( $#curr_data_name != 1 ) then
                echo "ERROR: No data file or multiple data files with ${data_stem} for bold run number ${bold_num}. Please check file existance or the bold list."
                exit 1;
            endif

            set data_names = "${data_names} ${curr_data_name}"

        end
        echo ${data_names} >> $output
        popd > /dev/null
    endif

    if( $folder == "surf" ) then
        pushd ${folder} > /dev/null

        foreach hemi (lh rh)
            set data_names = ""
            foreach bold_num ($bold)
                set curr_data_name = "${hemi}.${sub}_bld${bold_num}${data_stem}"
                set curr_data_name = "${sub_dir}/${sub}/${folder}/${curr_data_name}"
                if( $#curr_data_name != 1 ) then
                    echo "ERROR: No ${hemi} data file or multiple ${hemi} data files with ${data_stem} for bold run number ${bold_num}. Please check file existance or the bold list."
                    exit 1;
                endif

                set data_names = "${data_names} ${curr_data_name}"

            end
            if( $hemi == "lh" ) then
                echo ${data_names} >> ${lh_output}
            endif
            if( $hemi == "rh" ) then
                echo ${data_names} >> ${rh_output}
            endif
        end

        popd > /dev/null
    endif

    popd > /dev/null
end


popd > /dev/null
exit 0;


#############################
# parse arguments
#############################
parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
    set flag = $argv[1]; shift;

    switch($flag)
        # subjects directory
        case "-sd":
            if( $#argv == 0 ) goto arg1err;
            set sub_dir = $argv[1]; shift;
            breaksw

        # subject list
        case "-sub_ls":
            if( $#argv == 0 ) goto arg1err;
            set sub_list = $argv[1]; shift;
            breaksw

        # relative data folder, e.g. bold, vol, surf
        case "-folder":
            if( $#argv == 0 ) goto arg1err;
            set folder = $argv[1]; shift;
            breaksw

        # data stem
        case "-data_stem":
            if( $#argv == 0 ) goto arg1err;
            set data_stem = $argv[1]; shift;
            breaksw

        # output directory
        case "-out_dir":
            if( $#argv == 0 ) goto arg1err;
            set out_dir = $argv[1]; shift;
            breaksw

        # output stem
        case "-out_stem":
            if( $#argv == 0 ) goto arg1err;
            set out_stem = $argv[1]; shift;
            breaksw

        # Assumption for preprocessing pipeline
        case "-preproc_opt":
            if( $#argv == 0 ) goto arg1err;
            set preproc_opt = $argv[1]; shift;
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

if( $sub_dir == "" ) then
    echo "ERROR: subjects directory not specified."
    exit 1;
endif

if( $sub_list == "" ) then
    echo "ERROR: subject list not specified."
    exit 1;
endif

if( $folder == "" ) then
    echo "ERROR: relative data folder not specified."
    exit 1;
endif

if( $data_stem == "" ) then
    echo "ERROR: data stem not specified."
    exit 1;
endif

if( $out_dir == "" ) then
    echo "output directory not specified."
    exit 1;
endif

if( $out_stem == "" ) then
    echo "output stem not specified."
    exit 1;
endif

goto check_params_return;


###########################
# Error message
###########################
arg1err:
    echo "ERROR: flag $flag requires one argument"
    exit 1;


###########################
# Usage exit
###########################
usage_exit:

    echo ""
    echo "USAGE: CBIG_fMRI_create_data_list.csh"
    echo ""
    echo "    -sd           sub_dir     : subjects directory"
    echo "    -sub_ls       subjects    : subjects list"
    echo "    -folder       folder      : relative data folder, e.g. bold, vold, surf..."
    echo "    -data_stem    data_stem   : data stem (after *bld002, with extension)"
    echo "    -out_dir      out_dir     : output directory"
    echo "    -out_stem     out_stem    : output list stem"
    echo "    -preproc_opt  preproc_opt : assumption of preprocessing approach, choose from 'old' and 'new', 'old' means procsurffast file structure, 'new' means CBIG_preproc_fMRI_preprocess file structure. Default is 'new'"
    echo ""

    if ( $PrintHelp == 0 ) exit 1
    echo $VERSION

    cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1

#-------- Everything below is printed as part of help --------#
BEGINHELP

  Given the subjects directory, the subjects list, a relative data folder (e.g. bold, vol, surf), and the data stem, this function create a list (or two lists for surface data) where each line is one subject with different runs.
  This function assumes the fMRI data are preprocessed by CBIG_preproc_fMRI_preprocess.csh.
  The user needs to specify the output directory and output list stem.

