#!/bin/bash
#

# CBIG_check_whether_function_used_in_other_functions.sh $input_function_name $folder
# Search for all instances of a function name inside a given folder and give a warning


input_function_name=$1
folder=$2

# function to join array into string
function join_str { local IFS="$1"; shift; echo "$*"; }

# function to get lines from a file into a comma separated string 
function get_lines_from_file() {
  files_arr=()
  while IFS= read -r line # Read a line
  do
    files_arr+=("$line")
  done < "$1"
  result=$(join_str , "${files_arr[@]}")
  echo $result
}

# define prefix and file extensions that will be ignored
EXCLUDED_EXTENSIONS=$(get_lines_from_file "$CBIG_CODE_DIR/setup/check_function_format/excluded_extensions.txt")

# in a given folder, find all lines containing a function name without $PREFIX (first `grep`) and ignore lines containing the function name with the $PREFIX already preprended (third `grep`)
all_matches=$(grep -IHnRl $input_function_name --exclude=\*.{$EXCLUDED_EXTENSIONS} --exclude="${input_function_name}" --exclude-from="$CBIG_CODE_DIR/setup/check_function_format/excluded_files.txt" $folder)

if [ ! -z "$all_matches" ]; then   
   echo "$all_matches"
fi