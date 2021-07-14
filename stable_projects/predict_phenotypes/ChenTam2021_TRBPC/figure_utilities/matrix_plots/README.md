# Matrix plots
This folder contains code to generate predictive feature matrices for the project Task-Rest Behavioral Prediction in Children

# Usage
* Step 1 `CBIG_TRBPC_aggregate_relevance.m` : This function will aggregate the relevance values (and average across folds/seeds & normalize) across the behaviors and fMRI conditions and generate `.mat` files containing average the relevance values (per score)
* Step 2 `CBIG_TRBPC_plot_avg_relevance.m` : This function will average the activations according to a user-specified clustering/grouping of the predicted behaviors and produce a `.mat` file of those average weights (per cluster) - requires `relevance_vectors.mat` from Step 1
* Step 3 `CBIG_TRBPC_plot_masked_permuted_matrix.m` : This function will mask out the non-significant edges/network blocks and plots the significant edges/network blocks per behavioral cluster and fMRI condition. Also generates a `matrices_cat_masked_sig.mat` containing these masked matrices
* Step 4 `CBIG_TRBPC_edges_conjunction.m` : This function will do a conjunction analysis and plot matrices of the conjunction and relevance maps after conjunction. also generates  `.csv` files containing average scores for each network block in a symmetric matrix for the chord diagrams. Requires `matrices_cat_masked_sig.mat` from Step 3
* Step 5 `CBIG_TRBPC_chord_matrix_merge_pos_neg.m` : This function will put the significant positive and significant negative blocks together into one figure for the relevance matrix and chord diagrams. Requires `mat_relevance_conjunction_masked.mat` from Step 4
* Step 6 `CBIG_TRBPC_subcortical_heatmap.m` : This function will create heatmaps for the subcortical ROIs - requires `vec_relevance_sum_conjunction_subcortical.mat` from Step 4
* Step 7 `CBIG_TRBPC_supp_relevance_ind.m` : This function will plot the predictive feature matrices for all individual predicted behaviors

# This folder also contains the following functions which are called by the functions in the previous section
* `CBIG_TRBPC_rearrange_matrix.m` : This function will simply rearrange a 419x419 matrix into major networks (no plotting) - This work is derived from Siyi Tang
* `CBIG_TRBPC_plot_fc_matrix.m` : This function simply plots a 419x419 matrix (that is already arranged into the major networks and subcortical regions) - This work is derived from Siyi Tang
* `CBIG_TRBPC_newmatrix_2_oldmatrix_order.m` : This function will transform a matrix that was rearranged in the order from `CBIG_rearrange_matrix.m` back to the original ordering (needed for surface projection) - This work is derived from Siyi Tang
* `CBIG_TRBPC_FC_vector_2_mat.m` : This function will transform a vector into a matrix - by Jianzhong Chen
* `CBIG_TRBPC_generate_block_pairs_18networks.m` : This function will output the names and indices of 17 networks + the subcortical (as in `CBIG_TRBPC_rearrange_matrix.m`) and all the pairs of networks
* `CBIG_TRBPC_avg_network_pairs_relevance.m` : This function will average values within network blocks for a 419x419 matrix that is already arranged according to the scale 18 networks (17 cortical + 1 subcortical) order



