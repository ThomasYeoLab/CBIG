This folder contains the scripts to run linear ridge regression to predict behavioral measures for each dataset. Note that this folder only contains the wrapper scripts specific for the GSP and HCP. General linear ridge regression scripts are stored in `$CBIG_CODE_DIR/utilities/matlab/predictive_models/LinearRidgeRegression` since it is widely used in CBIG lab. Details about how to use our general linear ridge regression package are [here](https://github.com/ThomasYeoLab/CBIG/blob/master/utilities/matlab/predictive_models/LinearRidgeRegression/README.md).

- `GSP/scripts/CBIG_LiGSR_LRR_workflowGSP.sh`: the top-level wrapper script for a particular behavior in the **GSP** dataset, calling `GSP/scripts/CBIG_LiGSR_LRR_workflowGSP.m`. Type `CBIG_LiGSR_LRR_workflowGSP.sh` in command line for the details of usage.

- `GSP/scripts/CBIG_LiGSR_LRR_workflowGSP.m`: the matlab function that integrates the whole workflow of linear ridge regression for a particular behavior in the **GSP** dataset, calling general functions in `$CBIG_CODE_DIR/utilities/matlab/predictive_models/LinearRidgeRegression`.

- `HCP/scripts/CBIG_LiGSR_LRR_workflowHCP.sh`: the top-level wrapper script for a particular behavior in the **HCP** dataset, calling `HCP/scripts/CBIG_LiGSR_LRR_workflowHCP.m`. Type `CBIG_LiGSR_LRR_workflowHCP.sh` in command line for the details of usage.

- `HCP/scripts/CBIG_LiGSR_LRR_workflowHCP.m`: the matlab function that integrates the whole workflow of linear ridge regression for a particular behavior in the **HCP** dataset, calling general functions in `$CBIG_CODE_DIR/utilities/matlab/predictive_models/LinearRidgeRegression`.
