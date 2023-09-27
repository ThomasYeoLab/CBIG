#!/bin/bash
#
# CBIG_check_whether_function_used_in_other_functions_wrapper.sh $file_path "silent"
# Wrapper function to search for all instances of a function name 
# inside a predefined set of directories 
# Input can be the function name, a file name, or full path of a file


file_name=$(basename "$1")
verbose=$2

# if the input function is matlab, the function name has no extension
# if the input function is bash or cshell, the function name has extension
if [[ $file_name == *.m ]]; then
    function_name="${file_name%.*}"
elif [[ $file_name == *.sh ]] || [[ $file_name == *.csh ]]; then
    function_name="$file_name"
fi

# search for all instances of a function name
DIRECTORY_NAMES=("utilities" \
"stable_projects" \
"data/templates")
all_matches=""
i=0
for name in "${DIRECTORY_NAMES[@]}"
do
    i=$((i+1))
    directory="$CBIG_CODE_DIR/$name"
    if [[ "$verbose" != "silent" ]]; then
    echo "* Looking into $directory"
    fi
    matches=$($CBIG_CODE_DIR/setup/check_function_format/CBIG_check_whether_function_used_in_other_functions.sh $function_name $directory)
    if [ ! -z "$matches" ]; then
        if [ $i == 1 ]; then
            all_matches="$matches"
        else
            all_matches="$all_matches $matches"
        fi
    fi
done
if [ ! -z "$all_matches" ]; then
   echo "[WARNING] You have changed $file_name which has been used in other functions."
   echo "$all_matches"
   echo "We advise you to inform and double check with CBIG lab memebers who are using above functions."
fi
