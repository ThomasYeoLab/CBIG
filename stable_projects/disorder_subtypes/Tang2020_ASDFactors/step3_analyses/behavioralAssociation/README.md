This folder performs canonical correlation analysis (CCA) between factor loadings and behavioral scores to investigate the associations between latent factors and behavioral symptoms.

Figure 5 was created by functions in this folder.

----
## What Does Each File Do?
1. `CBIG_ASDf_CCA_genInputs.m` prepares inputs for CCA. E.g., constructs exchangeability block for within-site permutation, constructs regressors etc.
2. `CBIG_ASDf_CCA_factorBehavior.m` performs CCA with permutation test (within-site permutation) between factor loadings and behavioral scores.
3. `CBIG_ASDf_CCA_plotBar.m` creates bar plots for behavioral scores structure coefficients (i.e., Pearson's correlation between behavior CCA loading & original behavioral scores).
4. `CBIG_ASDf_CCA_plotScatter.m` creates scatter plots for factor CCA loading vs behavior CCA loading.
5. `CBIG_ASDf_CCA_wrapper.m` is the wrapper function to perform CCA between each factor loading & each group of behavioral scores (see Methods section in our paper for details) for K = 3.