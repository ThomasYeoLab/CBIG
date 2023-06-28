# Parametric Mean Field Model Analysis for Schaefer 100 parcellation
* `input` this folder contains the data to perform the estimation process of pMFM on Schaefer 100 parcellation. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the Schaefer 100 results of the paper. 
    * `CBIG_pMFM_basic_functions_Schaefer100.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_step1_training_Schaefer100.py` this file is the wrapper to perform parameter estimation process. The process performs 10 random initializations and each random initialization takes 500 iterations. In total, there are 5000 candidate parameter sets.
    * `CBIG_pMFM_step2_validation_Schaefer100.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_Schaefer100.py`. The 5000 candidate parameter sets are evaluated in the validation set to obtain top 10 candidate parameter sets.
    * `CBIG_pMFM_step3_test_Schaefer100.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_Schaefer100.py`. The script generates 1000 simulations for each parameter set and computes the averaged cost.
    * `CBIG_pMFM_step4_generate_simulated_fc_fcd.py` this file is the wrapper to generate simulated FC and FCD based on MFM and estimated model parameters from `CBIG_pMFM_step3_test_Schaefer100.py`
    * `CBIG_pMFM_step5_generate_STDFCD_correlation_Schaefer100.m` this file is the function to compute the empirical and simulated SWSTD-FCD correlation
    * `CBIG_pMFM_step6_SWSTD_state_Schaefer100.m` this file is the function to perform SWSTD state analysis.
    * `CBIG_pMFM_step7_perturbation_analysis.py` this file is the wrapper to perform perturbation analysis
    * `CBIG_pMFM_step8_gene_expression_analysis_desikan.m` this file is the wrapper to perform gene expression analysis 

# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Schaefer100_parcellation/scripts` and follow the step order from step 1 to step 7 to run the scripts.
* For Python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts. 
* For Matlab scripts, directly run the scripts with proper inputs.