# Append `CBIG_` to functions
In matlab or other programming languages, it is easy to have two functions which have same names. For example, matlab built-in function `sort` and your own function `sort`. In order to avoid this conflict between different functions from different sources, our lab, CBIG (Computational Brain Imaging Group), append `CBIG_` to all functions written by our group. If you do not append `CBIG_` to your functions, you can use the scripts in this folder to append `CBIG_` automatically.

## Steps 
These scripts will   
1) Check whether a given function has prefix `CBIG_`. If this function does not have `CBIG_`, the scripts will append `CBIG_` to the function name.  
2) Loop all files in folders "utilities" "stable_projects" "setup" "data/templates" to find scripts which use the given function.   
PS: The searching will exclude files with extension in file `excluded_extensions.txt` and exclude files in file `excluded_files.txt`.  
3) Give a prompt to ask you whether you want to replace the given function in those scripts with new name `CBIG_XXX`. If you want to replace, click `y`. 

## Usage
* Once you `git commit` in your `CBIG_private` repo, the git hook will automatically check all functions in your commit and prepend `CBIG_` to all those functions without `CBIG_`.
* If you want to use these functions by yourself, you can just type 
```
CBIG_prepend_prefix_to_function_name.sh $function_name
```


