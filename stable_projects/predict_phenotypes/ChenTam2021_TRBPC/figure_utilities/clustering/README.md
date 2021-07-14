# Clustering
This folder contains code to generate the similarity matrices among predctive network features for the project Task-Rest Behavioral Prediction in Children.
More specifically, this folder contains Jupyter notebooks written in Python.

# What does this folder contain?
* `CBIG_TRBPC_python_env.txt` : A specification file to build a conda environment that is identical to the one that was used to produce the original figures
* `multikernel_weight_corr_clustering_datadriven.ipynb` : this notebook will perform a hierarchical clustering of the relevance maps/predictive network feature matrices (in vector form). Requires output from `../matrix_plots/CBIG_TRBPC_aggregate_relevance.m`
* `multikernel_weight_corr_clustering_hyp.ipynb` : this notebook will group the relevance maps/predictive network feature matrices (in vector form) from `plot_weight_multikernel.m` based on hypothesis-driven clusters. Requires output from `../matrix_plots/CBIG_TRBPC_aggregate_relevance.m`
* `similarity_fmri_clusters.ipynb` : this notebook will generate similarity matrices for the average clustered relevance maps across fMRI states (e.g. cognition rest, cognition SST, etc). Requires output from `../matrix_plots/CBIG_TRBPC_plot_avg_relevance.m`
* `similarity_behaviour.ipynb` : this notebook will perform a hierarchical clustering of the actual behavioral scores.

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
```bash
# Activate conda environment
source activate py_trbpc

# Open up a notebook
jupyter notebook notebooknamehere
```
