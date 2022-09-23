# HCP Behaviour Prediction
This folder contains scripts pertaining to behaviour prediction in the HCP.

## Data and splits
The features and behaviours used for regression models can be found in [HCP_repo/data](insert_github_link).
Similarly, the training and test splits can be found in [HCP_repo/splits](insert_github_link).

## Regression procedure
We perform a 10-fold cross validation procedure in the HCP, with 60 random initializations. 
Age and sex were regressed from each behaviour score. As age of the HCP subjects requires restricted access, it is not shared on our HCP repository. Please apply for permission [here](HCP_link) and access it from the HCP website. The list of subjects that we use in this study can be found in [HCP_repo/data/behaviours](insert_github_link). 
