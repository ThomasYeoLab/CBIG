# Control analysis by only using RSFC principal gradient for the parameters parametrization
* Normally, the model parameters of pMFM are parameterized by both RSFC gradient and T1w/T2w ratio. Here in this control analysis, we only consider the RSFC gradient.
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this will be used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_step1_training_gradient.py` this file is the wrapper to perform parameter estimation process. The process performs 10 random initializations and each random initialization takes 500 iterations. In total, there are 5000 candidate parameter sets.
    * `CBIG_pMFM_step2_validation_gradient.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_gradient.py`. The 5000 candidate parameter sets are evaluated in the validation set to obtain top 10 candidate parameter sets.
    * `CBIG_pMFM_step3_test_gradient.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_gradient.py`. The script generates 1000 simulations for each parameter set and computes the averaged cost.


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Gradient_only/scripts` and follow the step order from step 1 to step 3 to run the scripts.
* For python scripts, please first activate the environment (e.g. `source acitvate pMFM`) and then directly run the scripts.