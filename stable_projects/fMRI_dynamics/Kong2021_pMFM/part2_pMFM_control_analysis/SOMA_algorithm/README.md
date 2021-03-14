# Control analysis of using SOMA algorithm
* In the main results of the paper, we utilize CMA-ES to optimize model parameters. Here, the control analysis here utilizes SOMA and does not perform as well as CMA-ES
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_test_SOMA_training.py` this file is the wrapper estimate model parameters using SOMA algorithm


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/SOMA_algorithm/scripts` and run the script `CBIG_pMFM_test_different_window.py`
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.