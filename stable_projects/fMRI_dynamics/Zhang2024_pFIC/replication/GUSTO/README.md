## GUSTO
This folder contains the input data, configuration files and reference outputs for the analyses that utilized the GUSTO dataset.

---
### Overview
The GUSTO dataset is mainly used for the purpose of showing the replicability and generalizability of the PNC cognitive effect results. The procedure was largely the same as that of the PNC cognitive effect analysis. Here, we considered 5 cognitive scores and derived a composite score by taking the their first principal component. However, due to the much smaller sample size compared to PNC, we did not further divide the high- and low-performance groups into sub-groups. Thus, only 1 set of E/I ratio was estimated for the high- or low-performance group each. The statistical significance was established by randomly assigning the composite score to each participant and re-spliting them into high- and low-performance groups and estimating the E/I ratio. This random permutation procedure was done for 100 times, creating a null distribution. 

The results were consistent with the PNC cognitive effect analysis results. We observed that the E/I ratio was lower for the high-performance group than low-performance group. The regional-level effect sizes (Cohen's d) were larger in association regions compared with sensory regions.

---
### input

The `input` folder contains the following files that are used to train, validate the model and estimate the E/I ratio.

* FC and FCD of the training, validation and test set of low and high-performance group. Note that generating these inputs requires access to the participant's fMRI time series. (see `util/GUSTO/CBIG_generate_GUSTO_subject_level_FC_FCD.m` and `util/GUSTO/CBIG_generate_GUSTO_group_level_FC_FCD.m` for more details)
* SC of the training, validation and test set. Same as HCP's SCs.
* Subject lists of training, validation and test set showing subject ID for easy replication/reference with corresponding age and cognitive score.
* RSFC gradient and T1w/T2w ratio (`myelin`) are the same as those used for HCP training.

---
### config
This folder contains the configuration files used for estimating E/I ratio of each age group or low/high-performance groups.

---
### script
Wrapper functions to replicate the analyses can be found under this folder.

---
### reference_output

This folder includes the results used to generate Figure 6 of the manuscript. The final E/I ratio estimates for each age group or low/high-performance group are saved as `EI.csv` under respective `test` folder. Values represent E/I ratio estimate for each ROI according to the Desikan parcellation (See `util/general/Desikan_72.txt` for reference, in particular, ROI #1, #5, #37 and #41 are removed as they correspond to the medial wall).

---
### How can I replicate these results?

To exactly replicate the results requires consistent hardware and software versions. All the analyses were run on RTX3090 GPU with CUDA 11.7 and the python environment file is located here `replication/config/CBIG_pFIC_python_env.yml`. If the aim is to replicate the results, make sure that you are using the same GPU version and Python environment prior to running the replication wrapper functions.
 
As an example, run the wrapper function as follows to replicate the estimated E/I ratio
```
bash CBIG_pFIC_wrapper_GUSTO.sh
```
Typical run-time using one RTX3090 GPU is about 10 hours. (This script does not include running the permutation test)

