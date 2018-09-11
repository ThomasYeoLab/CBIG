## Examples

Here, we provide a example to make sure that users outside CBIG lab can run the code in `Wang2018_MFMem/step2_simulation/` correctly.
In this example, we provide our empirical functional connectivity (FC) and structure connectivity (SC) matrix and estimated parameters which are in Desikan parcellation. These FC and SC are stored in `FCSC_Desikan68_Raphael_Wang.mat` under `data/` folder. The estimated parameters are stored in the `Example_Estimated_Parameter.mat` under `data/` folder.
The user can define the following variables as input for the script `CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_example.m`:

* `save_dir_input`: Output file directory. The users should define their own file saving directory.
* `FCSC_file_name`: The `.mat` file containing the testing FC and SC matrix.
* `parameter_file_name`: The `.mat` file containing the estimated parameters obtained from step 1. 

To run the our example, the user should run the following commands:

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Wang2018_MFMem/examples/step2_simulation_example/
```
Start Matlab, in Matlab command window, run the following commands 
```
CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_example(save_dir_input,'FCSC_Desikan68_Raphael_Wang','Example_Estimated_Parameter')
```

If the user would like to apply our code on their own data, please note that our scripts only applicable for FC and SC in Desikan parcellation. Note that Desikan parcellation has 36 ROIs in each hemisphere, we ignore the data in ROIs `unkown` and `corpuscallosum`. Therefore, in our case, there are 34 ROIs in each hemisphere.

Step1: Make sure all your input empirical FC and SC matrix are in Desikan parcellation and stored in a `.mat` file. This `.mat` file shoud contain two variables FC matrix called `FC_test` and SC matrix called `SC_test`.

Step2: Make sure your input estimated parameters are obtained in step 1 and stored in a `.mat` file. This `.mat` file shoud contain a variable called `Para_E`.

Step3: The above `.mat` file should be moved into the `Wang2018_MFMem/examples/step2_simulation_example/data/` folder.

Step4: The user can define the following variables as input for the script `CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_example.m`:

* `save_dir_input`: Output file directory. The users should define their own file saving directory.
* `FCSC_file_name`: The `.mat` file containing the testing FC and SC matrix.
* `parameter_file_name`: The `.mat` file containing the estimated parameters obtained from step 1. 

----

## References
Wang P, Kong R, Kong XL, Liegeois R, Orban C, Deco G, Van den Heuvel M, Yeo BT. Inversion of a Large-Scale Circuit Model Reveals 
a Cortical Hierarchy in the Dynamic Resting Human Brain

----

## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).
