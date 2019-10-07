This folder performs bootstrapping, re-estimates factors from the bootstrapped samples, and computes bootstrapped z-scores to obtain statistically significant RSFC patterns. Furthermore, this folder also includes functions to plot sinificant RSFC patterns (Figure 2B), RSFC patterns shared across all factors (Figure 4B), as well as RSFC patterns unique to each factor (Figure S2).

Figures 2B-2C, 4 and S2 are relevant to this folder.

----
## What Does Each File Do?
1. `CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly.m` plots 18x18 matrix by averaging/summing within/between networks in the 419x419 matrix.
2. `CBIG_ASDf_Plot400Schaefer19Subcor17Networks_thresholded.m` plots thresholded 419x419 matrix. Values below the threshold are set to 0.
3. `CBIG_ASDf_getBootstrappedSamples.m` performs bootstrapping and writes the corresponding RSFC into documents for factor estimate.
4. `CBIG_ASDf_bootstrappedEst_wrapper.sh` re-estimates factors for the bootstrapped samples.
5. `CBIG_ASDf_computeBootstrapZScores_wrapper.m` computes bootstrapped z-scores based on the re-estimated factors.
6. `CBIG_ASDf_plotFactorsThresholded_wrapper.m` plots the thresholded RSFC patterns (Figure 2B) based on the bootstrapped z-scores.
7. `CBIG_ASDf_plotConjunctionUniqueMaps_wrapper.m` finds RSFC patterns shared across factors as well as unique to each factor, and plots them (Figure 4B, S2).
