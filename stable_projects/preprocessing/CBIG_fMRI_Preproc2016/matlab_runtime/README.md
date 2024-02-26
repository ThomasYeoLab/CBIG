**Folder Structure:**
This folder contains 2 folders:

`examples` folder: contains example input txt file containing names of functions-to-be-compiled.

`scripts` folder: contains `CBIG_preproc_matlab_runtime_compile_utilities.m` script, which uses MATLAB mcc function to compile functions.
For more details on the script, please refer to the description within the script itself.

Running `CBIG_preproc_matlab_runtime_compile_utilities.m` generates a 3rd folder `utilities` which contains the executables (for example, `CBIG_preproc_FCmetrics`) and a mcc log file (`CBIG_matlab_runtime_compile_utilities.log`).


**How to compile MATLAB functions:**
1) Prepare a txt file `utilities2compile.txt` whereby each line corresponds to the name of each to-be-compiled-function.
For example, you intend to compile CBIG_fMRI_Preproc2016 utility functions `CBIG_preproc_FCmetrics.m`, and `CBIG_preproc_create_ROI_regressors.m`.
Then `utilities2compile.txt` should have two lines: "CBIG_preproc_FCmetrics.m" followed by "CBIG_preproc_create_ROI_regressors.m". (without quotes)
You can refer to `${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/matlab_runtime/utilities/examples/input/utilities2compile.txt`
for a concrete example.

2) Run `${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/matlab_runtime/CBIG_preproc_matlab_runtime_compile_utilities.m`, with `utilities2compile.txt` as an input to the `util_file` argument. Refer to `CBIG_preproc_matlab_runtime_compile_utilities.m` for more details. Note that the licensed MATLAB version is required to run `CBIG_preproc_matlab_runtime_compile_utilities.m` for compilation.


**How to run MATLAB Runtime version of CBIG_fMRI_Preproc2016:**
1) Copy and paste the following line into your `.bashrc` file (should be under path `~/.bashrc`). The line is `export MATLAB_RUNTIME_DIR=<path_to_MATLAB_Runtime>/MCR_R2018b/v95`, which creates an environment variable pointing to MATLAB Runtime. Note that after editing your `.bashrc`, you need to relaunch your terminal for this environment variable to be created.

2) Append `-matlab_runtime` to $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh command. For example, `$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh -s 123 -output_d <path_to_output_dir> -anat_s FreeSurfer -anat_d <path_to_reconall_output> -fmrinii <path_to_fmrinii>/fmrinii_123.txt -config <path_to_config>/config_123.txt -matlab_runtime`. Similarly, if you are running a single preprocessing step, the procedure is identical. (i.e. append `-matlab_runtime` to the end of command)
