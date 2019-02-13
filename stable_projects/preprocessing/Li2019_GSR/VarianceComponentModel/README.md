This folder contains codes to run variance component model on each dataset separately. 

Note that in our paper, a jackknife procedure was utilized to quantify the uncertainty of the variance estimates. For each behavioral phenotype, half the participants were randomly removed. The variance component model was fitted for the remaining participants. This jackknife procedure was repeated 1000 times, resulting in 1000 jackknife estimates of the explained behavioral variance.

GSP functions
==============
- `scripts/CBIG_LiGSR_LME_workflowGSP.sh`: the top-level wrapper script for the **GSP** dataset, calling `scripts/CBIG_LiGSR_LME_workflowGSP.m`. Type `CBIG_LiGSR_LME_workflowGSP.sh` in command line for the details of usage.

- `scripts/CBIG_LiGSR_LME_workflowGSP.m`: the matlab function that integrates computing functional similarity matrix (FSM), generating jackknife samples, and performing variance component model on each jackknife sample for the **GSP** dataset.

- `scripts/CBIG_LiGSR_explained_variance_GSP.m`: given a subject list and behavior list, this function performs variance component model specifically in the **GSP** dataset.

HCP functions
==============
- `scripts/CBIG_LiGSR_LME_workflowHCP.sh`: the top-level wrapper script for the **HCP** dataset, calling `scripts/CBIG_LiGSR_LME_workflowHCP.m`. Type `CBIG_LiGSR_LME_workflowHCP.sh` in command line for the details of usage.

- `scripts/CBIG_LiGSR_LME_workflowHCP.m`: the matlab function that integrates computing FSM, generating jackknife samples, and performing variance component model on each jackknife sample for the **HCP** dataset.

- `scripts/CBIG_LiGSR_explained_variance_HCP.m`: given a subject list and behavior list, this function performs variance component model specifically in the **HCP** dataset.

Function shared across datasets
==============

- `utilities/CBIG_LiGSR_compute_FSM_from_FC.m`: given a RSFC file, this function computes the FSM matrix based on the lower-triangular, off-diagonal entries.

- `utilities/CBIG_LiGSR_NchooseD_families.m`: given a subject list, this function generate the list of subjects to be removed for each jackknife sample, considering family structure.

- `CBIG_LiGSR_LME_cmp2pipe_allstats.m`: this function calculates the statistics mentioned in Li et al. (under review) related to variance component model.

- `CBIG_LiGSR_del_d_jack_cmp2pipelines.m`: it computes 
  1. the percentage improvement of the explained variance (baseline+GSR pipeline versus baseline pipeline)
  2. the jackknife mean and variance of improved explained variance in (i). 
  
  This function is called by `CBIG_LiGSR_LME_cmp2pipe_allstats.m`.

- `CBIG_LiGSR_PosNeg_jackIQR_cmp2pipe.m`: it computes (1) the number of behavioral measures with entire interquartile range (IQR) of explained behavioral variance **difference** above 0 (or below 0); (2) the number of behavioral measures with the median of explained behavioral variance **difference** above 0 (or below 0).
