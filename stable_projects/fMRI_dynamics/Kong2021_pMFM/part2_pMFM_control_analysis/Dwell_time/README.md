# Control analysis of computing Dwell time of empirical and simulated FCD mean
* `scripts` this folder contains the scripts to generate the dwell time results of the paper. 
    * `CBIG_pMFM_count_func.m` this function is used to compute the dwell time if you give a FCD mean
    * `CBIG_pMFM_dwell_time_empirical.m` this file is the wrapper to compute the dwell time of empirical FCD mean
    * `CBIG_pMFM_dwell_time_simulated.m` this file is the wrapper to compute the dwell time of simulated FCD mean


# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Dwell_time/scripts` and run the script `CBIG_pMFM_dwell_time_empirical.m` or `CBIG_pMFM_dwell_time_simulated.m`
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.