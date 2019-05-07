This folder contains code for 10 fold cross validation. 1. We match the factor orders across folds. 2. We stack factor loading of 10 test foldes together and correlate atrophy loadings with behavior loadings.   

Figure 6 is created by functions in this folder.

----

## What Does Each File Do?
1. `CBIG_MMLDA_10foldCV_wrapper.m` is the wrapper function to do 10 fold cross validation. 
2. `CBIG_MMLDA_hungarian_match_factors.m` do hungarian match between factors from differnet folds.