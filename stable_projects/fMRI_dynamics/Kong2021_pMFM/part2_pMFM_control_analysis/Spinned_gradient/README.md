# Control analysis of computing results based on model parameters parameterized by spinned FC gradient and T1w/T2w
* `scripts` this folder contains the scripts to generate the results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_step1_training_SpGrad_SpT1T2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level spinned first FC gradient and spinned T1w/T2w
    * `CBIG_pMFM_step2_validation_SpGrad_SpT1T2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_SpGrad_SpT1T2.py`.
    * `CBIG_pMFM_step3_test_SpGrad_SpT1T2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_SpGrad_SpT1T2.py`.
    * `CBIG_pMFM_step4_training_SpGrad.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level spinned first FC gradient and original T1w/T2w
    * `CBIG_pMFM_step5_validation_SpGrad.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step4_training_SpGrad.py`.
    * `CBIG_pMFM_step6_test_SpGrad.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step5_validation_SpGrad.py`. 
    * `CBIG_pMFM_step7_training_SpT1T2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level original first FC gradient and spinned T1w/T2w
    * `CBIG_pMFM_step8_validation_SpT1T2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step7_training_SpT1T2.py`.
    * `CBIG_pMFM_step9_test_SpT1T2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step8_validation_SpT1T2.py`. 



# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Spinned_gradient/scripts` and run the scripts.
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.
