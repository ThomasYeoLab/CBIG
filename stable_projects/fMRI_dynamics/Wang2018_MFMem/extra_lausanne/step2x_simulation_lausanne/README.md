### Usage

This folder contains all matlab data, functions and scripts to use dynamic mean-field model to simulate the functional connectivity (FC) from estimated model parameters.
All matlab functions and scripts should be standalone but you need install matlab statistical toolbox. 

---

### Content
 * data/
 This sub-folder contains all necessary data (e.g. FC & SC).

 * lib/
 This sub-folder contains all functions.

   CBIG_mfm_rfMRI_mfm_ode1.m
   CBIG_mfm_rfMRI_nsolver_eul_sto_resLH.m
   CBIG_mfm_rfMRI_rfMRI_BW_ode1.m
   CBIG_mfm_rfMRI_simBOLD_downsampling.m

 * save/
 This sub-folder is reserrved to save the result.
 
 *  CBIG_mfm_rfMRI_run_1000_simulation_testGrp_resHigh.m 
   This function is to run 1000 times simluation of FC, useing estimated parameter and
   SC of test grp, simulation running in high time resolution t = 0.0001s.
    
 *  CBIG_mfm_rfMRI_run_1000_simulation_testGrp_resLow.m 
   This function is to run 1000 times simluation of FC, useing estimated parameter and
   SC of test grp, simulation running in high time resolution t = 0.01s.

 *  CBIG_mfm_rfMRI_run_1000_simulation_trainingGrp_resHigh.m 
   This function is to run 1000 times simluation of FC, useing estimated parameter and
   SC of training grp, simulation running in high time resolution t = 0.0001s.

 *  CBIG_mfm_rfMRI_run_1000_simulation_trainingGrp_resLow.m 
   This function is to run 1000 times simluation of FC, useing estimated parameter and
   SC of training grp, simulation running in high time resolution t = 0.01s.

For the content details of each sub-folder, please refer to the "readme.md" in each sub-folder seperately or read the comment inside each function & script. 

---

### Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).
