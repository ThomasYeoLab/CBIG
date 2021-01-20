### Usage

This folder contains all matlab data, functions and scripts to use dynamic mean-field model to simulate the functional connectivity (FC) obtained from resting-state fMRI. 
Inverse parameter estimation method is applied to fit the model parameter to the empiricial data. All matlab functions and scripts should be standalone but you need install matlab statistical toolbox. 

---

### Content
 * cluster/
 This sub-folder is reserved for the parallel computating. 
 
 * data/
 This sub-folder contains all necessary data (e.g. FC & SC).                
   
 * lib/
 This sub-folder contains all functions.
   
   CBIG_mfm_rfMRI_diff_P1.mat                             
   CBIG_mfm_rfMRI_diff_PC1.mat                
   CBIG_mfm_rfMRI_matrixQ.mat
   CBIG_mfm_rfMRI_mfm_ode1.mat
   CBIG_mfm_rfMRI_nsolver_eul_sto.mat
   CBIG_mfm_rfMRI_rfMRI_BW_ode1.mat
   CBIG_mfm_rfMRI_simBOLD_downsampling.mat
   CBIG_mfm_rfMRI_Trace_AXB.mat
   
 * save/
 This sub-folder is reserrved to save the result.
   
  
 *  CBIG_mfm_rfMRI_estimation_main.m
 This is the execution script, run it!


For the content details of each sub-folder, please refer to the "README.md" in each sub-folder seperately or read the comment inside each function & script.
----

### Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).

