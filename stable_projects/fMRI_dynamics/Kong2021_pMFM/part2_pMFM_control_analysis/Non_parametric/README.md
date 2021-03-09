# Control analysis without using model parameters parameterization
* Recurrent connection strength w, external input current I, and noise amplitude are allowed to be spatially heterogeneous across brain regions, but not constrained by T1w/T2w or FC gradient
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_step1_training_nonpara.py` this file is the wrapper to perform parameter estimation process. The process performs 10 random initializations and each random initialization takes 500 iterations. In total, there are 5000 candidate parameter sets.
    * `CBIG_pMFM_step2_validation_nonpara.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_nonpara.py`. The 5000 candidate parameter sets are evaluated in the validation set to obtain top 10 candidate parameter sets.
    * `CBIG_pMFM_step3_test_nonpara.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_nonpara.py`. The script generates 1000 simulations for each parameter set and computes the averaged cost.


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Non_parametric/scripts` and follow the step order from step 1 to step 3 to run the scripts.
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.