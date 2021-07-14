# Wrappers
This folder contains scripts/functions to generate figures for the project Task-Rest Behavioral Prediction in Children.
More specifically, these functions will generate all the figures that were written in MATLAB or R (e.g. predictive feature matrices, chord diagrams).
To generate the green-brown similarity matrices, which were written in Python, please see the notebooks in `../clustering`

## Usage
* `CBIG_TRBPC_pfm_surf_chord_wrapper.sh` : This will generate all the predictive feature matrices, cortical surface projections, subcortical heatmaps, and chord diagrams.
* `CBIG_TRBPC_matrix_plots_wrapper.m` : This MATLAB function utilizes functions in `../matrix_plots` and generates all the matrix plots of the predictive feature matrices. The outputs of this function are needed by other functions in `../chord`, `../surface_plots`, `../clustering`.
* `CBIG_TRBPC_replicate_matrix_plots.m`: This function calls `CBIG_TRBPC_matrix_plots_wrapper.m` to replicate replicate all the matrix in Chen & Tam 2021 paper. This function is called by the wrapper `CBIG_TRBPC_pfm_surf_chord_wrapper.sh`