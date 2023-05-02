# Interpretation codes
This folder contains all used to interpret feature importance from single-type-feature prediction models in the ABCD and HCP through Haufe's inversion. Scripts specific to each dataset are available in their respective folders.

## Folder structure in `HCP` and `ABCD` folders
The interpretation codes for each dataset are arranged as such:
1. `CBIG_MMP_<dataset>_interpretation_wrapper`: Wrapper script to run feature importance analysis for all models.
2. `CBIG_MMP_<dataset>_Haufe.sh`: Shell script for job to be submitted to the scheduler.
3. `CBIG_MMP_<dataset>_Haufe.m`: Matlab script calculating the Haufe-inversion for single-type-feature prediction models.
