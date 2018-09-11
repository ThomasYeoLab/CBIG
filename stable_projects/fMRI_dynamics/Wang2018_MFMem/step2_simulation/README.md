
## Content

* Sub-folder `data/` contains all necessary data (e.g. FC & SC).

  ```
  Desikan34ROIsperHemisphere_regionDescription.mat
  Example_Estimated_Parameter.mat  
  FCSC_Desikan68_Raphael_Wang.mat                                             
  ```
 
* `CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resHigh.m` 
This function will 1000 simluations of FC with estimated parameter and
SC of test group, simulation will run in high time resolution t = 0.0001s.
    
* `CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow.m` 
This function will 1000 simluations of FC with estimated parameter and
SC of test group, simulation will run in low time resolution t = 0.01s.

* `CBIG_MFMem_rfMRI_run_1000_simulation_trainingGrp_resHigh.m` 
This function will 1000 simluations of FC with estimated parameter and
SC of training group, simulation will run in high time resolution t = 0.0001s.

* `CBIG_MFMem_rfMRI_run_1000_simulation_trainingGrp_resLow.m` 
This function will 1000 simluations of FC with estimated parameter and
SC of training group, simulation will run in low time resolution t = 0.01s.

* In the `Wang2018_MFMem/lib/` contains all functions.

  ```
  CBIG_MFMem_rfMRI_mfm_ode1.m
  CBIG_MFMem_rfMRI_nsolver_eul_sto_resLH.m
  CBIG_MFMem_rfMRI_rfMRI_BW_ode1.m
  CBIG_MFMem_rfMRI_simBOLD_downsampling.m
  ```

----

## Usage

To use this function, the user should define the following variable as input for the 4 scripts.

* `save_dir_input`: Output file directory. The users should define their own file saving directory.

To run our function, the user should run the following commands:

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Wang2018_MFMem/step2_simulation/
```
Start Matlab, in Matlab command window, run the following commands (take one of the four functions as an example)
```
CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow(save_dir_input)
```

----

## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).
