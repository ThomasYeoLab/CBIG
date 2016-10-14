#!/bin/bash
#

# CBIG_prepend_prefix_to_function_name.sh $input_function_name $folder
# Search for all instances of a function name (without prefix) inside a given folder
# and replace them by the new function name with the prefix prepended, if prompted

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
PREFIX="CBIG_"
EXCLUDED_EXTENSIONS=$(get_lines_from_file "$CBIG_CODE_DIR/setup/check_function_format/excluded_extensions.txt")

# define old and new function name
if [[ $input_function_name == *"$PREFIX"* ]]; then
  old_function_name=${input_function_name#$PREFIX}
  new_function_name=$input_function_name
else
  old_function_name=$input_function_name
  new_function_name="${PREFIX}$input_function_name"
fi

# in a given folder, find all lines containing a function name without $PREFIX (first `grep`) and ignore lines containing the function name with the $PREFIX already preprended (third `grep`)
all_matches=(`grep -IHnR $old_function_name --exclude=\*.{$EXCLUDED_EXTENSIONS} --exclude-from="$CBIG_CODE_DIR/setup/check_function_format/excluded_files.txt" $folder | cut -d: -f1,2`)
last_file=""
for match in "${all_matches[@]}"
do
  # present the match
  file=$(echo $match | cut -d: -f1)
  line_number=$(echo $match | cut -d: -f2)
  line=$(sed "${line_number}q;d" $file)
  if [[ $line  != *$new_function_name* ]]; then
    echo ""
    old_line="  > Current line $line_number : $line"
    proposed_line="  > Proposed line $line_number: $line"
    proposed_line=${proposed_line/$old_function_name/$new_function_name}

    if [ "$file" != "$last_file" ]; then
      echo ""
      echo "Found this file: $file"
    fi

    # replace instances of the old function name with the new function name (having $PREFIX prepended)
    export GREP_COLOR='1;37;41'
    echo "$old_line" | grep --color='auto' -E "$old_function_name|$"
    export GREP_COLOR='1;32'
    echo "$proposed_line" | grep --color='auto' -E "$new_function_name|$"
    read -r -p "  Replace? [y/N] " response </dev/tty
    if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
    then
      replacement_line=${line/$old_function_name/$new_function_name}
      sed -i "${line_number}s~.*~$replacement_line~g" $file
      new_line=$(sed "${line_number}q;d" $file)
      new_line="  >> NEW line $line_number:     $new_line"
      echo "$new_line" | grep --color='auto' -E "$new_function_name|$"
    else
      echo "  >> Skipped"
    fi

    last_file=$file
  fi
done
