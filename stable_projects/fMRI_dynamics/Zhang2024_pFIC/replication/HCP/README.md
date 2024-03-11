## HCP
This folder contains the input data, configuration files and reference outputs for the analyses that utilized the HCP dataset.

---
### Overview
We used the HCP dataset to evaluate the efficacy of **p**arameterized **F**eedback **I**nhibition **C**ontrol (**pFIC**) model. We found that the model can fit to the empirical FC and FCD better (i.e., lower cost in the test set) when the model parameters (i.e., `wEE`, `wEI` and `sigma`) are parameterized by **both** T1w/T2w ratio (`myelin`) and the first principal RSFC gradient (`rsfc_gradient`). Thus, this folder contains the results of comparing 4 different ways to parameterize the model parameters. Illustrating with `wEE`

Method 1: pFIC: `wEE` = `a`  + `b` x `myelin` + `c` x `rsfc_gradient`\
Method 2: RSFC gradient only: `wEE` = `d`  +  `e` x `rsfc_gradient`\
Method 3: myelin only: `wEE` = `f`  +  `g` x `myelin`\
Method 4: homogeneous: `wEE` = `h`

where `a` - `h` are unknown constants that are optimized by CMA-ES during the model training process.

---
### input

The `input` folder contains the following files that are used to train, validate and test the model or to replicate certain figures shown in the manuscript.

* FC and FCD of the training, validation and test set (see `util/HCP/CBIG_generate_HCP_group_level_FC_FCD.m` for more details of how they are generated).
* SC of the training, validation and test set (see `util/general/CBIG_generate_group_level_SC.m` for more details).
* myelin and RSFC gradient used for parameterization, `vector_of_zeros.csv` is used as a placeholder to replace myelin and/or RSFC gradient (see Method 2 - 4 above).
* `TC.mat` is the parcellated time course from subject 103111 session 1-LR, which is used to generate Fig 2E of the manuscript.

---
### config
This folder contains the configuration files used for running the analyses indicated by their filenames.

---
### script
The CMA-ES initial parameters for evaluating Method 2, 3 and 4 are slightly different from those for Method 1. Instead of setting every CMA-ES parameter as a variable, which makes the scripts unnecessarily general and overly complicated, specific scripts were written for certain specific analyses. 

Wrapper functions to replicate the analyses can be found under this folder.

---
### reference_output

This folder includes the results used to generate Figure 2G of the manuscript.

---
### How can I replicate these results?

To exactly replicate the results requires consistent hardware and software versions. All the analyses were run on RTX3090 GPU with CUDA 11.7 and the python environment file is located here `replication/config/CBIG_pFIC_python_env.yml`. If the aim is to replicate the results, make sure that you are using the same GPU version and Python environment prior to running the replication wrapper functions.
 
As an example, run the wrapper function as follows to replicate the pFIC results.
```
cd script
bash CBIG_pFIC_wrapper_HCP.sh
```
Typical run-time using a RTX3090 GPU is about 25 hours.

