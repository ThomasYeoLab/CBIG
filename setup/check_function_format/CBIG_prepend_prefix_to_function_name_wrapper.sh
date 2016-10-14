#!/bin/bash
#

# CBIG_prepend_prefix_to_function_name_wrapper.sh $file_path
# Wrapper function to search for all instances of a function name (without prefix)
# inside a predefined set of directories and replace them by the new function name
# with the prefix prepended, if prompted
# Input can be the function name, a file name, or full path of a file

file_name=$(basename "$1")
function_name="${file_name%.*}"
verbose=$2

DIRECTORY_NAMES=("utilities" "stable_projects" "setup" "data/templates")
for name in "${DIRECTORY_NAMES[@]}"
do
  directory="$CBIG_CODE_DIR/$name"
  if [[ "$verbose" != "silent" ]]; then
    echo "* Looking into $directory"
  fi
  $CBIG_CODE_DIR/setup/check_function_format/CBIG_prepend_prefix_to_function_name.sh $function_name $directory
done
