#!/bin/bash
#

# CBIG_check_format_and_license_in_all_functions.sh
# Wrapper function to check format and license of all functions inside a predefined set of directories

DIRECTORY_NAMES=("utilities" "stable_projects" "setup" "data")
EXTENSIONS_TO_CHECK=("m" "sh" "csh" "c" "cpp" "py" "r")

for name in "${DIRECTORY_NAMES[@]}"
do
  directory="$CBIG_CODE_DIR/$name"
  for extension in "${EXTENSIONS_TO_CHECK[@]}"
  do
    command="find $directory -name '*.$extension'"
    eval $command | while read file_path; do
      echo "> Checking $file_path"
      # TO-DO: check license for source files other than Matlab files
      if [[ "$extension" == "m" ]]; then
        $CBIG_CODE_DIR/setup/check_license/CBIG_check_license_matlab_file.sh $file_path
      fi
      $CBIG_CODE_DIR/setup/check_function_format/CBIG_prepend_prefix_to_function_name_wrapper.sh $file_path silent
    done
  done
done
