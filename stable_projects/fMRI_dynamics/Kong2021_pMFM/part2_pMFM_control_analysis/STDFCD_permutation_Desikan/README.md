# Control analysis of performing permutation test on SWSTD-FCD correlation
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_step1_generate_permutation_order_desikan.m` this is the function that generate 10000 permutation order for Desikan results.
    * `CBIG_pMFM_step2_STDFCD_permutation_correlation_desikan.m` this is the function that compute 10000 SWSTD-FCD correlation based on the permutation order generated in step 1. The statistical test is also performed based on the null distrtibution .


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/High_resolution/scripts` and follow the step order to run the scripts.
* For Matlab scripts, you can run it directly in terminal.