# Clustering
This folder contains code to generate the similarity matrices among predctive network features for the project Task-Rest Behavioral Prediction in Children.
More specifically, this folder contains python scripts to generate the figures.

# What does this folder contain?
* `CBIG_TRBPC_python_env.txt` : A specification file to build a conda environment that is identical to the one that was used to produce the original figures
* `multikernel_PFM_similarity_plots.py` : this function can perform a hierarchical clustering of the relevance maps/predictive network feature matrices (in vector form). Requires output from `../matrix_plots/CBIG_TRBPC_aggregate_relevance.m`
* `similarity_fmri_clusters.py` : this function will generate similarity matrices for the average clustered relevance maps across fMRI states (e.g. cognition rest, cognition SST, etc). Requires output from `../matrix_plots/CBIG_TRBPC_plot_avg_relevance.m`
* `similarity_behavior.py` : this function will perform a hierarchical clustering of the actual behavioral scores.

# Installation
Because the code is written in Python, you need to install Python and the required packages first.

If you are using our `circuv`, I suggest the following installation method.
1. Install `miniconda`
To install `miniconda`, you can check `setup/python_env_setup`
2. Create and build a conda environment with `CBIG_TRBPC_python_env.txt`

```bash
# Create and build a new conda environment called py_trbpc
conda create --name py_trbpc --file CBIG_TRBPC_python_env.txt
```

# Usage
After the environment successfully created, check `../wrappers/CBIG_TRBPC_pfm_surf_chord_wrapper.sh` for example usage.