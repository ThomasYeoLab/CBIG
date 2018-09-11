
## Content

* Sub-folder `data/` contains all necessary data (e.g. FC & SC).

   ```
  Desikan34ROIsperHemisphere_regionDescription.mat
  FCSC_Desikan68_Raphael_Wang.mat                               
  saved_original_random_generator.mat
   ```

* `CBIG_MFMem_rfMRI_estimation_main_only4p.m`
This is the execution script for estimating the 4 parameters in the model, make sure that you are under `../Wang2018_MFMem/extra_4p/step1_estimation_4p` to run the function.

* All functions needed by the two scripts are listed as follow. The functions can be found in the folder `../Wang2018_MFMem/lib`

    ```
   CBIG_MFMem_rfMRI_diff_P1.mat                             
   CBIG_MFMem_rfMRI_diff_PC1.mat                
   CBIG_MFMem_rfMRI_matrixQ.mat
   CBIG_MFMem_rfMRI_mfm_ode1_4p.mat
   CBIG_MFMem_rfMRI_nsolver_eul_sto_4p.mat
   CBIG_MFMem_rfMRI_rfMRI_BW_ode1.mat
   CBIG_MFMem_rfMRI_simBOLD_downsampling.mat
   CBIG_MFMem_rfMRI_Trace_AXB.mat
    ```

----

## Usage

* This folder contains all matlab data and scripts to use dynamic mean-field model with only 4 parameters to simulate the functional connectivity (FC) obtained from resting-state fMRI. 

* Inverse parameter estimation method is applied to fit the model parameter to the empiricial data. All matlab functions and scripts should be standalone but you need install matlab statistical toolbox. 

* To run the functions `CBIG_MFMem_rfMRI_estimation_main_only4.m`, the user should define the following variables in the scripts:
1. `N_core`: Number of CPU cores for calculation, the default value is 30
2. `EstimationMaxStep`: Number of estimation setp, the default value is 500
3. `save_dir_input`: Output file directory. The user should define their own file-saving directory.

To run the our example, the user should run the following commands:

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Wang2018_MFMem/extra_4p/step1_estimation_4p/
```
Start Matlab, in Matlab command window, run the following commands 
```
CBIG_MFMem_rfMRI_estimation_main_only4p
```

----

## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).

