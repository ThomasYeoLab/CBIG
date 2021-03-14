# Control analysis of using high time resolution during fMRI simulation
* Normally, the time resolution for fMRI simulation in pMFM is 0.01s. Here in this control analysis, we change the time resolution to 0.001s.
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_test_high_resolution.py` this file is the wrapper to conduct test process based on 0.001s time resolution.


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/High_resolution/scripts` and run the script `CBIG_pMFM_test_high_resolution.py`
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.