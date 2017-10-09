#!/bin/bash
#

# CBIG_replace_old_with_new_function_name.sh $old_function_name $new_function_name
# Wrapper function to search for all instances of old function name
# inside a predefined set of directories and replace them by the new function name,
# if prompted
# Input can be the function name, a file name, or full path of a file

old_function_name=$(basename "$1")
old_function_name="${old_function_name%.*}"
new_function_name=$(basename "$2")
new_function_name="${new_function_name%.*}"


DIRECTORY_NAMES=("utilities" "stable_projects" "setup" "data/templates")
for name in "${DIRECTORY_NAMES[@]}"
do
  directory="$CBIG_CODE_DIR/$name"
  if [[ "$verbose" != "silent" ]]; then
    echo "* Looking into $directory"
  fi
  $CBIG_CODE_DIR/setup/replace_old_with_new_func_name/CBIG_replace_old_with_new_function_name.sh $old_function_name $new_function_name $directory
done
