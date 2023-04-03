# Comparison between gradients and parcellations for RSFC behavioral prediction

References
==========
+ Kong et al., [Comparison Between Gradients and Parcellations for Functional Connectivity Prediction of Behavior](https://doi.org/10.1016/j.neuroimage.2023.120044). NeuroImage. 2023

----

Background
====

Resting-state functional connectivity (RSFC) is widely used to predict behavioral measures. To predict behavioral measures, representing RSFC with parcellations and gradients are the two most popular approaches. Here, we compare parcellation and gradient approaches for RSFC-based prediction of a broad range of behavioral measures in the Human Connectome Project (HCP) and Adolescent Brain Cognitive Development (ABCD) datasets. Two regression algorithms were utilized for prediction: linear ridge regression (LRR) and kernel ridge regression (KRR).

<img width="427" alt="image" src="https://user-images.githubusercontent.com/20438248/228835479-02918113-8526-428d-ada6-9244704fb001.png">

----

Prediction Models
====

In this work, we used two prediction models LRR and KRR.
+ KRR scripts can be found here: https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/KernelRidgeRegression

+ LRR scripts can be found here: https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/LRR_fracridge

To see how to use the scripts, please check the `example` folder. To explore the experiment setup pf the current study, check the `Usage` section below.   

----

Data Release
====

Data from the following list will be released. However, as these files require a large amount of space, we are still exploring the best way to release them. Please check back for updates.

+ Schaefer2018 parcellations (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 ROIs) and corresponding RSFC matrices.
+ Kong2021 parcellations (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 ROIs) and corresponding RSFC matrices.
+ sICA maps (50, 100, 200, 300 ICA components) and corresponding RSFC matrices.
+ LocalGrad (Laumann et al., 2015; single gradient map)
+ PrincipalGrad (Margulies et al., 2016; 1, 5, 10, 20, 40, 80, 100 gradients)

The following data are available:

**Schaefer2018** 

+ Parcellations with different resolutions: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations

**Kong2021** 
+ Individual-level parcellations and FC (Pearson corr) with different resolutions for HCP dataset: https://github.com/ThomasYeoLab/Kong2022_ArealMSHBM
+ Individual-level parcellations with different resolutions for ABCD dataset: https://dx.doi.org/10.15154/1528046

**sICA**
+ ICA spatial maps with different resolutions for HCP data has been released in HCP website: http://humanconnectome.org/storage/app/media/documentation/s1200/HCP1200-
DenseConnectome+PTN+Appendix-July2017.pdf
+ sICA maps with different resolutions for ABCD dataset can be found here: https://dx.doi.org/10.15154/1528046

----

Code Release
====

The KRR and LRR prediction models were released in another directory:

https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models

Here we released the scripts to optimize resolutions:
+ `CBIG_GradPar_LRR_frac_optimize_res_wrapper.m`
+ `CBIG_GradPar_KRR_optimize_res_wrapper.m` 

**Examples**

We provide detailed examples of above two steps in **`examples`** folder. **We highly recommended the users to go through the example tutorial first**.

**Download**

To download the version of the code that is last tested, you can either

- visit this link: [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.28.0-Kong2023_GradPar](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.28.0-Kong2023_GradPar)

or

- run the following command, if you have Git installed

```
git checkout -b Kong2023_GradPar v0.28.0-Kong2023_GradPar
```

----

Usage
====

The code utilized in this study include two main steps:

+ Generate KRR and LRR prediction results for each parcellation/gradient approach and each resolution.

+ Once the prediction results are generated, we can use the `CBIG_GradPar_LRR_frac_optimize_res_wrapper.m` and `CBIG_GradPar_KRR_optimize_res_wrapper.m` script to find the optimal resolution for each parcellation/gradient approach.

**Step1: Perform prediction using LRR and KRR for each resolution**

The LRR_fracridge and KRR scripts can be found here:

+ LRR_fracridge: https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/LRR_fracridge

+ KRR: https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/KernelRidgeRegression

To check how to use the scripts, please check the `README.md` files in the above two links.

**Step2: Find optimal resolution for prediction**

To optimize resolution for each behavioral measure, the preditction results from step1 are needed. The prediction settings have to be the same across all resolutions. These parameters were required for the wrapper script. Check the description of the wrapper script for how to construct the input parameters.

+ LRR_fracridge: `CBIG_GradPar_LRR_frac_optimize_res_wrapper.m`

+ KRR: `CBIG_GradPar_KRR_optimize_res_wrapper.m`


----

Updates
=======
- Release v0.28.0 (03/04/2023) Initial release of Kong2023_GradPar.

----

Bugs and Questions
====
Please contact Ru(by) Kong at roo.cone@gmail.com.
