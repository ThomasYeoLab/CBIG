## Alprazolam
This folder contains the input data, configuration files and reference outputs for the analyses that utilized the Alprazolam dataset.

---
### Overview
The Alprazolam dataset is used to validate our pFIC approach, which is used to estimate whole-cortex excitation/inhibition (E/I) ratio. 

In this dataset, participants underwent a double-blind, placebo-controlled study using the benzodiazepine alprazolam. Each participant completed two identical experimental sessions approximately 1 week apart. In one session, participants were given a 1-mg dose of alprazolam, and in the other, they were given an identical appearing placebo. Pharmacologically, **alprazolam** is a GABA-agonist that **increases inhibitory activities**. This means that the E/I ratio estimated from the alprazolam session should be lower than that from the placebo session. 

Our results show that, using pFIC model, estimated E/I ratio is statistically significantly higher in alprazolam session than placebo session. The statistical significance is established through a permutation test. Also, the estimated regional E/I ratio correlates with regional benzodiazepine receptor (BZR) density measured using PET. 

---
### input

The `input` folder contains the following files that are used to train, validate and test the model or to replicate certain figures shown in the manuscript.

* FC and FCD of the training, validation and test set of both drug and placebo sessions. Note that the generation of these inputs require access to the participant's fMRI time series. (see `util/Alprazolam/CBIG_generate_Alprazolam_subject_level_FC_FCD.m` and `util/Alprazolam/CBIG_generate_Alprazolam_group_level_FC_FCD.m` for more details)
* SC of the training, validation and test set. They are the same as HCP's SCs. 
* RSFC gradient and T1w/T2w ratio used during training for paramterization. 
* Subject lists of training, validation and test set showing subject ID for easy replication/reference.
* `roi_list_0.5.txt` is a binary vector indicating which ROIs (based on Desikan-Killiany parcellation) are retained for analysis. This is because the original fMRI acquisition has a limited FOV (See Figure S2), thus only regions with >50% ROI coverage are included (indicated by `1`). [**Note**: SC, RSFC gradient and T1w/T2w ratio, except during extrapolation, only have 42 ROIs. That is after masking out the regions with <50% ROI coverage]
* `BZR_density` is a 68-by-1 vector. Each entry represents regional benzodiazepine receptor density measured using PET (Norgaard et al., 2021).

---
### config
This folder contains the configuration files used for estimating E/I ratio of drug session and placebo session. A configuration file template used for permutation test is also provided.

---
### script
Wrapper functions to replicate the analyses can be found under this folder.

---
### reference_output

This folder includes the results used to generate Figure 3 of the manuscript. Additionally, a pre-generated permutation-testing-based E/I ratio null distributions for drug and placebo sessions are included.

---
### How can I replicate these results?

To exactly replicate the results requires consistent hardware and software versions.  All the analyses were run on RTX3090 GPU with CUDA 11.7 and the python environment file is located here `replication/config/CBIG_pFIC_python_env.yml`. If the aim is to replicate the results, make sure that you are using the same GPU version and Python environment prior to running the replication wrapper functions.
 
As an example, run the wrapper function as follows to replicate the estimated E/I ratio of the drug session.
```
cd script
bash CBIG_pFIC_wrapper_drug.sh
```
Typical run-time using a RTX3090 GPU is about 12 hours.

