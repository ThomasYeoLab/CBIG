# Add license in Matlab functions
We use MIT license for all codes in our lab, CBIG (Computational Brain Imaging Group). This script is used to automatically format your Matlab functions and add the MIT license.
PS: The scripts can only add license to Matlab functions for now. If you want to add MIT license in other bash or python scripts, you have to add it manually.

## Steps 
These scripts will   
1) Move the comment block below the function name  
2) Check and add function claim in comment block  
3) Check and add the MIT license below comment block

## Usage
* Once you `git commit` in your `CBIG_private` repo, the git hook will automatically check all functions in your commit and you can choose whether you want to move the comment block or add the license. 
* If you want to use these functions by yourself, you can just type 
```
CBIG_check_license_one_folder.sh $directory
```


