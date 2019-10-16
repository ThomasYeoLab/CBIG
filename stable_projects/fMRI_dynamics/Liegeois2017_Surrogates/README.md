## Reference

Raphael Liegeois, Timothy Laumann, Avi Snyder, Juan Zhou and B. T. Thomas Yeo. [**Interpreting Temporal Fluctuations In Resting-State Functional Connectivity MRI**](https://www.sciencedirect.com/science/article/abs/pii/S1053811917307516). NeuroImage, 163:437-55, 2017.

----
## Background

We show that (1) the space of stationary models is bigger than expected, and includes hidden Markov models (HMM), (2) Phase & autoregressive randomizations test for stationarity, linearity, Gaussianity, (3) Resting-state fMRI is mostly stationary, linear, and Gaussian. (4) 1st order autoregressive (AR) model encodes static & one-lag functional connectivity, and (5) 1st order AR model explains sliding window correlations very well, better than HMM.

----

## Code Release

This folder contains the code used to generate the results presented in our paper.

CBIG_RL2017_Figure4_toy_ex_HMM.m           reproduces the results presented in Figure 4 of the referenced paper
CBIG_RL2017_Figure7_toy_ex_bi_vs_multi.m   reproduces the results presented in Figure 7 of the referenced paper

Other results of the paper were obtained by combining autoregressive and phase randomization procedures. These procedures are implemented in CBIG_RL2017_get_AR_surrogate.m and CBIG_RL2017_get_PR_surrogate.m. Different options are detailed in the headers of the corresponding codes. 

To download the exact code utilized for the paper, you can either

- visit this link:
[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.15.3-Update_proj_refs_and_add_KRR_LITE](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.15.3-Update_proj_refs_and_add_KRR_LITE)

or

- run the following command, if you have Git installed
 
```
git checkout -b Liegeois2017_Surrogates v0.15.3-Update_proj_refs_and_add_KRR_LITE
```

----

## Updates

- Release v0.4.0 (17/09/2017): Initial release of Liegeois2017_Surrogates.
- Release v0.4.1 (17/09/2017): Update this README.md.
- Release v0.4.2 (09/10/2017): Add `CBIG_tested_config.sh` for Liegeois2017_Surrogates.
- Release v0.4.10 (01/02/2018):
    
    1. Add project-specific prefix `RL2017` for all scripts names.
    
    2. Add `config`, `unit_tests`, and `examples` folders.

- Release v0.15.3 (16/10/2019): Update README files.

----

## Bugs and Questions

Please contact Raphael Liegeois at Raphael.Liegeois@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
