This folder contains scripts about how to replicate results of this project. Notice that all filenames and directories in these scripts only work for CBIG lab.

# DATA
Input data used for replication is stored in `$CBIG_REPDATA_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM`
# REPLICATE ALL ANALYSIS
`./CBIG_pMFM_replication_all_wrapper.sh`
* This function will run all the analysis done in the paper by calling all other scripts in this folder that only replicate part of our analysis.
* To skip some part of the analysis, just comment out the lines in this function.
* To find the results of the replication, just go to `./partX_pMFM_XXX/ANALTSIS_NAME/results`
* To check if the replication is correct, compare with the results in `$CBIG_REPDATA_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/results`

# REPLICATE PART OF THE ANALYSIS
## Replicate the part 1 pMFM main analysis
`./CBIG_pMFM_replication_part1_pMFM_main.sh`
* This function will run the pMFM main analysis of our paper
* To find the results of the replication, just go to `./part1_pMFM_main/results`

## Replicate the part 2 pMFM control analysis
`./CBIG_pMFM_replication_part2_pMFM_control_analysis.sh`
* This function will run the pMFM control analysis of our paper
* To find the results of the replication, just go to `./part2_pMFM_control_analysis/ANALTSIS_NAME/results`
