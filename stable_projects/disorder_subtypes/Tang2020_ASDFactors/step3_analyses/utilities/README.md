This folder contains some helper functions that are called by other functions in multiple folders.

----
## What Does Each File Do?
1. `CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input.m` is a modified version of `$CBIG_CODE_DIR/utilities/matlab/figure_utilities/CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid.m`. It takes in a 419 x 419 matrix as input and plots it. This function is mainly used to visualize latent factors (i.e., E(RSFC patterns|Factor)).
2. `CBIG_ASDf_getSubData.m` is a helper function to extract participants' demographics/characteristics data from a spreadsheet.
2. `CBIG_ASDf_genRegressors.m` is a helper function to construct regressors (e.g., age, sex, head motion & sites)
3. `CBIG_ASDf_indivCorr2avgCorr.m` is a function to compute average correlation of participants' RSFC matrices from individual RSFC matrices.
