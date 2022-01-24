# Task-Rest Behavioral Prediction in Children (TRBPC)
# REFERENCE
* Chen, J., Tam, A., Kebets, V., Orban, C., Ooi, L.Q.R., Marek, S., Dosenbach, N., Eickhoff, S., Bzdok, D., Holmes, A.J. and Yeo, B.T., 2020. Shared and unique brain network features predict cognition, personality and mental health in childhood. BioRxiv.

# BACKGROUND
We used resting-state and task functional connectivity to predict a wide range of behaviors, spanning cognition, mental health symptoms, and personality, in the Adolescent Brain Cognitive Development study. We found that combining resting-state and task functional connectivity improves the prediction of cognition and personality. We also found that cognition, mental health, and personality were predicted by distinct patterns of brain organization.


![network features](readme_figures/network_features.png)
# CODE RELEASE
## Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: https://github.com/ThomasYeoLab/Standalone_ChenTam2022_TRBPC

## Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

To download the version of the code that was last tested, you can either

visit this link: https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.21.4-Update_ChenTam2022_TRBPC

run the following command, if you have Git installed
git checkout -b ChenTam2022_TRBPC v0.21.4-Update_ChenTam2022_TRBPC

# USAGE
## Setup
1. Make sure you have installed: Matlab 2018b
2. Follow $CBIG_CODE_DIR/setup/README.md to setup the config file. Instead of using CBIG/setup/CBIG_sample_config.sh, you need to use $CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/replication/config/CBIG_TRBPC_tested_config.sh.

## Replication
* `replication` this folder contains scripts to replicate all the analysis in this paper. This is mainly for internal use as it uses data on our server that are not able to made public.

## Examples and unit tests
* `examples` this folder provides a toy example for our kernel regression and linear regression workflow, as well as the computation of predictive feature matrix
* `unit_tests` this folder runs codes in `examples` and check with the reference output

## Regression models
* `KRR_LpOCV` (Kernel Ridge Regression Leave-p-Out Cross Validation): this folder contains leave-p-out cross validation workflow for both multi-kernel ridge regression and single-kernel ridge regression.
* `LRR_LpOCV` (Linear Ridge Regression Leave-p-Out Cross Validation): this folder contains leave-p-out cross validation workflow for linear ridge regression.
## Interpret the regression models
* This depends on the results from regression models
* `PFM` (Predictive-Feature Matrix): this folder contains functions to compute the predictive-feature matrix for both kernel ridge regression and linear ridge regression.
## Permutation test for significance
* This depends on the results from the regression models
* `permutation` this folder contains functions for the permutation test of kernel regression and predictive feature matrix
## Generate figures 
* `figure_utilities` this folder contains the scripts to generate the figures in our paper given you have all the results
## Post-hoc analysis
* `analysis` this folder contains functions to perform post-hoc analyses used in our paper



# UPDATES
* Release v0.21.4 (24/01/2022): Add codes for new analyses during the paper revision
* Release v0.21.0 (15/07/2021): Initial release of ChenTam2022_TRBPC project

# BUGS and QUESTIONS
Please contact Jianzhong Chen at chenjianzhong1996@gmail.com, Angela Tam at angela.tam08@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.