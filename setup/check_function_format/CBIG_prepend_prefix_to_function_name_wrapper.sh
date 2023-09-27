#!/bin/bash
#

# CBIG_prepend_prefix_to_function_name_wrapper.sh $file_path
# Wrapper function to search for all instances of a function name (without prefix)
# inside a predefined set of directories and replace them by the new function name
# with the prefix prepended, if prompted
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

DIRECTORY_NAMES=("utilities" \
"stable_projects/preprocessing" \
"stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering" \
"setup" \
"data/templates" \
"external_packages")
for name in "${DIRECTORY_NAMES[@]}"
do
    directory="$CBIG_CODE_DIR/$name"
    if [[ "$verbose" != "silent" ]]; then
        echo "* Looking into $directory"
    fi
    $CBIG_CODE_DIR/setup/check_function_format/CBIG_prepend_prefix_to_function_name.sh $function_name $directory
done
