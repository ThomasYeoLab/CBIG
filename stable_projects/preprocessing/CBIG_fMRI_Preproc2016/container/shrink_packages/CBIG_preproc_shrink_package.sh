#! /bin/sh

# This script is part of a process to shrink packages
# by removing redundant files not used in the preprocessing pipeline.
# In summary, the shrinking process involves:
# 1. running fMRI preprocessing unit tests and using strace to identify files accessed during execution.
# 2. we can moved some of the identified files to a temporary directory using this script
#    and test if the setup still works without them.
# 3. If the setup still works, we have successfully identified redundant files.
#    If the setup does not work, we can use this script move the files back
#    to their original location and try again with other identified files
# For more details, please refer to the README.md in the current directory.
#
# Specifically, this script takes two arguments: a relative path to a folder or file,
# and a flag ("shrink" or "restore").
# - If the flag is "shrink", it moves the specified folder or file from the APPS_DIR to TMP_APPS.
# - If the flag is "restore", it moves the specified folder or file back from TMP_APPS to APPS_DIR.
#
# Written by Tian Fang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

export CONTAINER_DIR=${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container
export APPS_DIR=${CONTAINER_DIR}/apps
export TMP_APPS=${CONTAINER_DIR}/tmp/apps
export LOG_FILE=${CONTAINER_DIR}/../unit_tests/single_subject/strace.log

# Function to handle the moving of files and directories
move_file_or_dir() {
    local src_dir=$1
    local dest_dir=$2
    local file_or_dir=$3

    if grep -q "$file_or_dir" "$LOG_FILE"; then
        if [ -d "$src_dir/$file_or_dir" ]; then
            for sub_item in $(ls "$src_dir/$file_or_dir"); do
                move_file_or_dir "$src_dir/$file_or_dir" "$dest_dir/$file_or_dir" "$sub_item"
            done
        fi
    else
        # make sure that the target file or directory does not exist in the $dest_dir
        if [ -e "$dest_dir/$file_or_dir" ]; then
            echo "Error: $dest_dir/$file_or_dir already exists"
            return
        fi
        mkdir -p "$(dirname $dest_dir/$file_or_dir)"
        mv "$src_dir/$file_or_dir" "$dest_dir/$file_or_dir"
        echo "Folder or file moved from $src_dir to $dest_dir"
    fi
}

# Check if the number of input arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <relative_path_to_folder_or_file> <restore_flag>"
    return
fi

# Check if the second input argument is a valid flag
if [ "$2" != "restore" ] && [ "$2" != "shrink" ]; then
    echo "Error: $2 is not a valid flag"
    return
fi

# Check if the second input argument is "shrink"
if [ "$2" == "shrink" ]; then
    move_file_or_dir "$APPS_DIR" "$TMP_APPS" "$1"
fi

# Check if the second input argument is "restore"
if [ "$2" == "restore" ]; then
    # Ensure the target folder does not exist
    if [ -d "$APPS_DIR/$1" ]; then
        echo "Error: $APPS_DIR/$1 already exists"
        return
    fi
    # Prepare the target folder in $APPS_DIR if it does not exist
    mkdir -p "$(dirname $APPS_DIR/$1)"
    # Move the folder or file back from $TMP_APPS to $APPS_DIR while preserving the relative file structure
    mv $TMP_APPS/$1 $APPS_DIR/$1
    echo "Folder or file moved back from $TMP_APPS to $APPS_DIR"
fi
