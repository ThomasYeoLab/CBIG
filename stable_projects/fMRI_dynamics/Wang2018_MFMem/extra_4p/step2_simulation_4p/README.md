
## Content

* Sub-folder `data/` contains all necessary data (e.g. FC & SC).

  ```
  Desikan34ROIsperHemisphere_regionDescription.mat
  Example_Estimated_Parameter.mat  
  FCSC_Desikan68_Raphael_Wang.mat                                             
  ```
    
* `CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_only4p.m` 
This function is to run 1000 times simluation of FC, using estimated parameter and
SC of test group, simulation running in high time resolution t = 0.01s.

* In the `Wang2018_MFMem/lib/` folder contains all functions related to this function.

  ```
  CBIG_MFMem_rfMRI_mfm_ode1_4p.m
  CBIG_MFMem_rfMRI_nsolver_eul_sto_resLH_4p.m
  CBIG_MFMem_rfMRI_rfMRI_BW_ode1.m
  CBIG_MFMem_rfMRI_simBOLD_downsampling.m
  ```

----

## Usage

* This folder contains all matlab data and scripts to use dynamic mean-field model to simulate the functional connectivity (FC) obtained from resting-state fMRI. 

* Inverse parameter estimation method is applied to fit the model parameter to the empiricial data. All matlab functions and scripts should be standalone but you need install matlab statistical toolbox. 

The user can define the following variable as input for the script `CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_only4p.m` :

* `save_dir_input`: Output file directory. The users should define their own file saving directory.

To run the our function, the user should run the following commands:

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Wang2018_MFMem/extra_4p/step2_simulation_4p/
```
Start Matlab, in Matlab command window, run the following commands (take one of the four functions as an example)
```
CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_only4p(save_dir_input)
```

----

## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).
