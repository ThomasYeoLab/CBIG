This folder generates Fig. S2.

This folder explores the possibility that two unknown factors in the (K+1)-factor model were subdivisions of an unknown factor in the K-factor model (whereas the other factors remained the same).

---

## What Does Each File Do?

1. `CBIG_wrapper.m` shows how to call the other functions.
2. `CBIG_factorHierachy.m` is the main function of this folder. It matches up factors from different models and computes correlations.
3. `CBIG_find_maxLikeRun.m` finds the random initialization with maximum likelihood for a given K.
4. `CBIG_effectiveNoSubjectsPerFactor.m` computes how many "effective subjects" we have for each factor. 
