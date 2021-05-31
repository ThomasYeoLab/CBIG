# Replication of Areal-leve Multi-session hierarchical Bayesian model (Areal MS-HBM)

The replication of MS-HBM contain several tasks:
1. The group priors estimation on 40 HCP data in fs_LR_32k and fsaverage6
2. The individual-level parcellation generation on two sets of data:
   - ABCD data in fsaverage6 surface space: 5 subjects, each subject has 4 sessions
   - HCP data in fs_LR_32k surface space: 755 subjects, each subject has 4 sessions
   
Notice that all filenames and directories in this replication **work for CBIG lab only**.

----

References
==========
+ Kong R, Yang Q, Gordon E, et al. [Individual-Specific Areal-Level Parcellations Improve Functional Connectivity Prediction of Behavior](https://doi.org/10.1093/cercor/bhab101). Cerebral Cortex. 2022

----

Data
====
In this replication, we will perform the algorithms on the functional connectivity profiles of 2 datasets:
+ **ABCD:** `$CBIG_REPDATA_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/ABCD_profile`
  
+ **HCP:** `$CBIG_ArealMSHBM_REP_HCP_DIR/??????/MNINonLinear/Results/rfMRI_REST?_??/postprocessing/MSM_reg_wbsgrayordinatecortex`

----

Run
====
In terminal, call the replication wrapper script:
```
$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/CBIG_ArealMSHBM_replication_wrapper.sh <output_dir>
```
In this script, it will submit several jobs:
1. estimating priors for fs_LR_32k using 3 HCP subjects with 4 sessions: 1 parent job + 3 children jobs
2. estimating priors for fsaverage6 using 3 HCP subjects with 4 sessions: 1 parent job + 3 children jobs
3. generate individual parcellation for 1 HCP subjects with 4 sessions: 1 job
4. generate individual parcellation for 1 ABCD subjects with 4 sessions: 1 job

Note that to fully replicate Kong2022, user can set the number of training subjects to be 40 instead of 3. We also only run two iterations to save time in this wrapper script, user can set the number of subjects and number of iterations as: all_sub=40; num_iter="" at line 27-28 in the wrapper script.Furthermore, we only generate individual parcellations for a single subject as examples.

### Group priors estimation

After running the above script, the estimated group priors will be saved into `priors` folder, the final results will be saved as `Params_Final.mat`, which contains a struct variable `Params`.

+ `<output_dir>/estimate_group_priors/HCP_fs_LR_32k/priors/gMSHBM/beta100`
+ `<output_dir>/estimate_group_priors/HCP_fsaverage6/priors/gMSHBM/beta5`

The estimated group priors include:
1) Inter-subject functional connectivity variability -- `Params.epsil`
2) Group-level connectivity profiles for each parcel -- `Params.mu`
3) Intra-subject functional connectivity variability -- `Params.sigma`
4) Spatial prior which denotes the probability of each parcel occurring at each location -- `Params.theta`

### Individual-level parcellations generation

The generated individual parcellations will be saved under:
+ `<output_dir>/generate_individual_parcellations/ABCD_fsaverage6/ind_parcellation_gMSHBM/test_set/2_sess/beta5`
+ `<output_dir>/generate_individual_parcellations/HCP_fs_LR_32k/ind_parcellation_gMSHBM/test_set/3_sess/beta50`


----

Check Results
====

### Group priors estimation
The user can compare 

+ `<output_dir>/estimate_group_priors/HCP_fs_LR_32k/priors/gMSHBM/beta100/Params_Final.mat`
+ `<output_dir>/estimate_group_priors/HCP_fsaverage6/priors/gMSHBM/beta5/Params_Final.mat`

with

+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/results/estimate_group_priors/HCP_fs_LR_32k/priors/gMSHBM/beta100/Params_Final.mat`
+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/results/estimate_group_priors/HCP_fsaverage6/priors/gMSHBM/beta5/Params_Final.mat`

### Individual-level parcellations generation
The user can compare

+ `<output_dir>/generate_individual_parcellations/HCP_fs_LR_32k/ind_parcellation_gMSHBM/test_set/3_sess/beta50/Ind_parcellation_MSHBM_sub1_w30_MRF30_beta50.mat`
+ `<output_dir>/generate_individual_parcellations/ABCD_fsaverage6/ind_parcellation_gMSHBM/test_set/2_sess/beta5/Ind_parcellation_MSHBM_sub1_w50_MRF50_beta5.mat`

with

+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/results/generate_individual_parcellations/HCP_fs_LR_32k/ind_parcellation_gMSHBM/test_set/3_sess/beta50/Ind_parcellation_MSHBM_sub1_w30_MRF30_beta50.mat`
+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication/results/generate_individual_parcellations/ABCD_fsaverage6/ind_parcellation_gMSHBM/test_set/2_sess/beta5/Ind_parcellation_MSHBM_sub1_w50_MRF50_beta5.mat`


----

Bugs and Questions
====
Please contact Ru(by) Kong at roo.cone@gmail.com.


