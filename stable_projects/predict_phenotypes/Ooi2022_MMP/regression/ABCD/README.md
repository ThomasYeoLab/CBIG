# HCP Behaviour Prediction
This folder contains scripts pertaining to behaviour prediction in the HCP.

## Data and splits
The features and behaviours used for regression models, as well as the training and test splits have been uploaded on the NDA website [NDA study link](NDA study link). Please apply for access to ABCD [here] on the NDA website to gain access to the study data.

## Regression procedure
We perform a 10-choose-3 site cross validation. Subjects from the ABCD are grouped into 10 sites, in each cross-validation fold we choose subjects from 3 sites as test set and the remaining subjects in the training set. This is repeated for all 120 possible combinations. 
Age and sex were regressed from each behaviour score.