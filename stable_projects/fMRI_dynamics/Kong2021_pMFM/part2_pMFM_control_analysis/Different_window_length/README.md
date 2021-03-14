# Control analysis of using different window length for sliding window
* Normally, the window length used to generate FCD is 83 TR. Here in this control analysis, we change the sliding window length to both 43 TR and 125 TR.
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the control analysis results of the paper. 
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_test_different_window.py` this file is the wrapper to conduct test process by using different sliding window length.


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Different_window_length/scripts` and run the script `CBIG_pMFM_test_different_window.py`
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.