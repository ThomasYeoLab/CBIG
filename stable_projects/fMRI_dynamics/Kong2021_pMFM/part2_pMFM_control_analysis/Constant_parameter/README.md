# Control analysis of setting model parameters (w, I, sigma) the same across regions
* Here, w refers to recurrent connections; I refers to external input; sigma refers noise amplitude
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_step1_training_conpara.py` this file is the wrapper to perform parameter estimation process. The process performs 10 random initializations and each random initialization takes 500 iterations. In total, there are 5000 candidate parameter sets.
    * `CBIG_pMFM_step2_validation_conpara.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_conpara.py`. The 5000 candidate parameter sets are evaluated in the validation set to obtain top 10 candidate parameter sets.
    * `CBIG_pMFM_step3_test_conpara.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_conpara.py`. The script generates 1000 simulations for each parameter set and computes the averaged cost.
    * `CBIG_pMFM_step4_generate_simulated_fc_fcd.py` this file is the wrapper to generate simulated FC and FCD based on MFM and estimated model parameters from `CBIG_pMFM_step3_test_conpara.py`
 

# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Constant_parameter/scripts` and follow the step order from step 1 to step 4 to run the scripts.
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.