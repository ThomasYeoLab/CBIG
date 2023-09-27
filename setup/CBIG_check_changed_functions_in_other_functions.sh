#!/bin/sh

###
# compare current branch with origin develop branch and list the changed files
###
curr_branch=$(git rev-parse --abbrev-ref HEAD)
files_to_be_checked=($(git diff --name-status $curr_branch..upstream/develop | awk '{ print $2 }'))

###
# define files to be checked
###
EXTENSIONS_TO_CHECK=("m" "sh" "csh")
DIRECTORIES_TO_CHECK=("stable_projects" "utilities" "external_packages")
EXCLUDED_FILES=("Surf2SurfGui.m" "Vol2SurfGui.m" "CBIG_tested_config.sh" "CBIG_tested_config.csh")

###
# check whether original function xxx of new function CBIG_xxx has been used in other functions
###
echo -e "\n==> Checking occurences of old function names (without CBIG prefix)"
for file_path in "${files_to_be_checked[@]}"
do
    file_name=( $(basename "$file_path") )

    # check whether file should be excluded
    file_in_excluded=0
    for excluded_file in "${EXCLUDED_FILES[@]}"
    do
        if [[ $file_name == $excluded_file ]]; then
            file_in_excluded=1
        fi
    done
    if [[ $file_in_excluded == 1 ]]; then
        break
    fi

    for ext in "${EXTENSIONS_TO_CHECK[@]}"
    do
        for directory in "${DIRECTORIES_TO_CHECK[@]}"
        do
            if [[ $file_path == $directory/* ]] && [[ $file_name == *.$ext ]]; then
                echo -e "\n  ==> Checking $file_path"
                $CBIG_CODE_DIR/setup/check_function_format/CBIG_prepend_prefix_to_function_name_wrapper.sh $file_path silent
            fi
        done
    done
done
echo "  [DONE]"

###
# check whether added/modifed/deleted functions have been used in other functions
###
echo -e "\n==> Checking occurences of added/modifed/deleted function names"
for file_path in "${files_to_be_checked[@]}"
do
    file_name=( $(basename "$file_path") )

    # check whether file should be excluded
    file_in_excluded=0
    for excluded_file in "${EXCLUDED_FILES[@]}"
    do
        if [[ $file_name == $excluded_file ]]; then
            file_in_excluded=1
        fi
    done
    if [[ $file_in_excluded == 1 ]]; then
        break
    fi

    for ext in "${EXTENSIONS_TO_CHECK[@]}"
    do
        for directory in "${DIRECTORIES_TO_CHECK[@]}"
        do
            if [[ $file_path == $directory/* ]] && [[ $file_name == *.$ext ]]; then
                echo -e "\n  ==> Checking $file_path"
                $CBIG_CODE_DIR/setup/check_function_format/CBIG_check_whether_function_used_in_other_functions_wrapper.sh $file_path silent
            fi
        done
    done
done
echo "  [DONE]"

###
# comments before git push
###
echo ""
echo "If you have renamed a function and want to change all instances of old function to the new function, you can use: "
echo "setup/replace_old_with_new_func_name/CBIG_replace_old_with_new_function_name_wrapper.sh"
echo ""
echo "If you have changed some files, please do the <git commit> again."

read -r -p "===> Do you still want to push current branch?[y/n]" response </dev/tty
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo "[PASSED]"
else
    echo "[FAILED] Abort pushing."
    exit 1
fi
