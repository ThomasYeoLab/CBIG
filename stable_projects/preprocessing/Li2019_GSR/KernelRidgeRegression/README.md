This folder contains the scripts to run kernel ridge regression to predict behavioral measures for each dataset. Note that this folder only contains the wrapper scripts specific for the GSP and HCP. General kernel ridge regression scripts are stored in `$CBIG_CODE_DIR/utilities/matlab/predictive_models/KernelRidgeRegression` since it is widely used in CBIG lab. Details about how to use our general kernel ridge regression package are [here](https://github.com/ThomasYeoLab/CBIG/blob/master/utilities/matlab/predictive_models/KernelRidgeRegression/README.md).

- `GSP/scripts/CBIG_LiGSR_KRR_workflowGSP.sh`: the top-level wrapper script for the **GSP** dataset, calling `GSP/scripts/CBIG_LiGSR_KRR_workflowGSP.m`. Type `CBIG_LiGSR_KRR_workflowGSP.sh` in command line for the details of usage.

- `GSP/scripts/CBIG_LiGSR_KRR_workflowGSP.m`: the matlab function that integrates the whole workflow of kernel ridge regression for the **GSP** dataset, calling general functions in `$CBIG_CODE_DIR/utilities/matlab/predictive_models/KernelRidgeRegression`.

- `HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh`: the top-level wrapper script for the **HCP** dataset, calling `HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.m`. Type `CBIG_LiGSR_KRR_workflowHCP.sh` in command line for the details of usage.

- `HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.m`: the matlab function that integrates the whole workflow of kernel ridge regression for the **HCP** dataset, calling general functions in `$CBIG_CODE_DIR/utilities/matlab/predictive_models/KernelRidgeRegression`.
