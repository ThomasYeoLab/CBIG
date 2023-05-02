# Run regression models

Functions in this folder provides functions to run the regression models. 

## Usage
We provide scripts for KRR, LRR, LASSO and random forest (RF) in their respective folders. 
*IMPORTANT:* Run the kernel ridge regression (KRR) first, as the other models will used the cross-validation subject groupings generated from the KRR.
Please read the instructions included in the scripts for directions on how to run the script. 

## Notes
* Please note that the LRR and LASSO depend on the same Elasticnet workflow as they are simply 2 cases of Elasticnet with alpha set to 0 and 1 respectively. 
* To run the RF model, please ensure you have the python dependencies installed. You can install `CBIG_ICCW_py_env.yml` under the `replication/config` folder.