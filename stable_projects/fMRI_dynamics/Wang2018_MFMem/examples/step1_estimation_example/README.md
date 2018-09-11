## Examples

Here, we provide a example to make sure that users outside CBIG lab can run the code in `Wang2018_MFMem/step1_estimation/` correctly.
In this example, we provide our empirical functional connectivity (FC) and structure connectivity (SC) matrix which are in Desikan parcellation. These FC and SC are stored in `FCSC_Desikan68_Raphael_Wang.mat` under `data/` folder.
The user can define the following variables inisde the script `CBIG_MFMem_rfMRI_estimation_main_example.m` (Line 79 to 82):

* `N_core`: Number of CPU cores for calculation. In this example, we set it to be 30. 
* `EstimationMaxStep`: Number of estimation setp. In this example, we set it to be 2.
* `save_dir_input`: Output file directory. The users should define their own output directory

To run the our example, the user should run the following commands:

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Wang2018_MFMem/examples/step1_estimation_example/
```
Start Matlab, in Matlab command window, run the following commands 
```
CBIG_MFMem_rfMRI_estimation_main_example
```

If the user would like to apply our code on their own data, please note that our scripts only applicable for FC and SC in Desikan parcellation. Note that Desikan parcellation has 36 ROIs in each hemisphere, we ignore the data in ROIs `unkown` and `corpuscallosum`. Therefore, in our case, there are 34 ROIs in each hemisphere.

Step1: Make sure all your input empirical FC and SC matrix are in Desikan parcellation (68 ROIs) and stored in a `.mat` file. This `.mat` file shoud contain two variables `FC_training` and `SC_training`.

Step2: The above `.mat` file should be moved into the `Wang2018_MFMem/examples/step1_estimation_example/data/` folder.

Step3: Modify Line 82 in script `CBIG_MFMem_rfMRI_estimation_main_example.m` to your `.mat` file:
```
datafile_name = <.mat file name containing FC and SC matrix>
```

Step4: Define the following variables in the script `CBIG_MFMem_rfMRI_estimation_main_example.m` (Line 79 to 81):

* `N_core`: Number of CPU cores for calculation.
* `EstimationMaxStep`: Number of estimation setp.
* `save_dir_input`: Output file directory. The users should define their own output directory

----

## References
Wang P, Kong R, Kong XL, Liegeois R, Orban C, Deco G, Van den Heuvel M, Yeo BT. Inversion of a Large-Scale Circuit Model Reveals 
a Cortical Hierarchy in the Dynamic Resting Human Brain

----


## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).

