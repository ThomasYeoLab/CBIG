# Random forest Leave-p-Out Cross-Validation workflow

Functions in this folder provides a wrapper function to perform leave-p-sites-out cross-validation workflow for RF. Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times. 

## Usage
To generate the RF prediction models, run the wrapper script `CBIG_ICCW_submit_RF_train.sh`. Please ensure that KRR was run already, as these scripts requires the param files generated before running KRR. The random forest scripts use python. Please ensure you have the python dependencies installed. You can install `CBIG_ICCW_py_env.yml` under the `replication/config` folder.
For the ease of reading the FC matrix, we have converted it to a csv. The user may wish to modify the script to read a mat file directly. The wrapper script runs the following:
1. `CBIG_ICCW_RF_train.py`: Python script that runs the random forest regression. 

To calculate the conditional variable importance of the RF models, run the wrapper script `CBIG_ICCW_submit_RF_cperm.sh`. Please ensure that the RF models have already finished running. 
For the ease of reading the FC matrix, we have converted it to a csv. The user may wish to modify the script to read a mat file directly. The wrapper script runs the following:
1. `CBIG_ICCW_RF_interpretation.py`: Python script that calculates the conditional variable importance

Please read the instructions included in the scripts for directions on how to run the script. 