# Replace old with new function name
If you want to rename a function, you may worry that other functions which call the old function may not work with new name. The scripts in this folder can check all functions and replace the function's old name with the new name automatically.

## Steps 
These scripts will   
1) Loop over all files under the folders "utilities" "stable_projects" "setup" "data/templates" to find scripts which use the old function.   
PS: The search will exclude files with extension listed under `excluded_extensions.txt` and exclude files listed under `excluded_files.txt`.  
2) Prompt to check if you want to replace the old function used in those scripts with the new function's name. Type `y` if you want to replace, `n` otherwise. 

## Usage
From a bash terminal, run:
```
CBIG_replace_old_with_new_function_name_wrapper.sh $old_function_name $new_function_name
```
