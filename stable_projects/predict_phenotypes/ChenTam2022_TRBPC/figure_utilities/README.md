# Figure utilities
This folder contains scripts and functions to generate figures for the project Task-Rest Behavioral Prediction in Children

## What does each folder do?
* `chord` : Contains code to generate chord diagrams. Written in R. Requires outputs from `matrix_plots/CBIG_TRBPC_chord_matrix_merge_pos_neg.m`
* `clustering` : Contains Jupyter notebooks to perform hierarchical clustering and create similarity matrices among the predictive network features. Written in Python. Requires outputs from `matrix_plots/CBIG_TRBPC_aggregate_relevance.m` and `matrix_plots/CBIG_TRBPC_aggregate_relevance.m` 
* `input` : Contains input files.
* `matrix_plots` : Contains code to generate the 419 x 419 predictive network feature matrices. Written in MATLAB. Requires output from `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/replication/CBIG_TRBPC_replication_all_wrapper.sh`
* `surface_plots` : Contains code to plot the relevance values on the cortical surface with the Schafer 400 parcellation. Requires output from `matrix_plots/CBIG_TRBPC_edges_conjunction.m`
* `wrappers` : Contains wrapper scripts to generate almost all the figures from the paper (predictive feature matrices, cortical surface projections, chord diagrams). Requires output from `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/replication/CBIG_TRBPC_replication_all_wrapper.sh`






