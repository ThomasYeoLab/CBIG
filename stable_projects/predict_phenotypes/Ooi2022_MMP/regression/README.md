# Regression codes
This folder contains all regression models used to predict behaviours in the ABCD and HCP. Scripts specific to each dataset are available in their respective folders, while scripts that both pipelines share are stored in the `utilities` folder.

## Folder structure in `HCP` and `ABCD` folders
The regression codes for each dataset are arranged as such:

### first_level
This folder contains scripts to run regression for a single feature type (e.g. resting state fmri only). The models include KRR, LRR and Elasticnet.
Each regression model has a .m file which you can run directly from matlab after parsing in the required arguments for each function. A .sh script is provided for each model to call the matlab script from your HPC if applicable.

### second_level
This folder contains scripts to run models that combine multiple feature types (e.g. combining resting and task fmri). The models multiKRR, and stacking. The stacking models require the KRR models from the first level to have been run already.
Each regression model .m file which you can run directly from matlab after parsing in the required arguments for each function. A .sh script is provided for each model to call the matlab script from your HPC if applicable.

### stats
This folder contains:
1. Scripts to create the null models for statistical tests, and calculate their p-value. To calculate the p-value, the original model must have been run already. 
2. Scripts to perform the pairwise comparison between models to see whether one performs better than another.

### utilities
This folder contains 3 main scripts. 
1. `CBIG_MMP_<dataset>_schedule_regression` \
This script initialises arguments for each regression type and stats scripts and submits them to the scheduler. 
2. `CBIG_MMP_<dataset>_read_model_results` \
This script navigates the ouput structure of each regression type and collates the final accuracy.
3. `CBIG_MMP_<dataset>_collate_result_wrapper` \
This script collates the final accuracy for all models run in the paper.

## Scripts in the `utilities` folder
This folder contains code used in both pipelines.
