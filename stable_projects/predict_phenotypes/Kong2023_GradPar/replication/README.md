# Replication of comparison between gradients and parcellations for RSFC behavioral prediction

The replication of the paper "Comparison Between Gradients and Parcellations for Functional Connectivity Prediction of Behavior" by Kong et al. (2023) contains the following steps:
1. Generate KRR prediction results for Kong2021 in HCP dataset with 100 and 200 resolutions. We will only run 3 splits.
2. Generate KRR prediction results for PrincipalGrad in ABCD dataset with 1 gradient. 
3. Generate LRR fracridge prediction results for Schaefer2018 in HCP dataset with 100 and 200 resolutions. We will only run 3 splits.
4. Generate LRR fracridge prediction results for sICA in ABCD dataset with 50 and 100 resolutions.
5. Optimize the resolution for Kong2021 KRR prediction in HCP dataset.
6. Optimize the resolution for sICA LRR fracridge prediction in ABCD dataset.

Please note that in the wrapper script `CBIG_GradPar_replication_wrapper.sh`, we only submit jobs for several representative prediction results in the paper to save time. To fully replicate the results in the paper, please modify the input parameters (e.g. project name and resolution) in the wrapper script:

+ Kong2021: 100,200,300,400,500,600,700,800,900,1000.
+ Schaefer2018: 100,200,300,400,500,600,700,800,900,1000.
+ PrincipalGrad: 1,5,10,20,40,60,80,100.
+ LocalGrad: 1
+ sICA: 50,100,200,300

For HCP dataset, please run 100 splits.

----

References
==========
+ Kong et al., [Comparison Between Gradients and Parcellations for Functional Connectivity Prediction of Behavior](https://doi.org/10.1016/j.neuroimage.2023.120044). NeuroImage. 2023

----

Data
====
The feature matrices of different approaches with different resolutions, behavior scores, covariates an other input parameters for HCP and ABCD datasets were saved in the CBIG server: 
`$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar`

----

Run
====
In terminal, call the replication wrapper script:
```
$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/replication/CBIG_GradPar_replication_wrapper.sh <output_dir>
```

The output structure:
+ `<output_dir>/KRR/HCP/Kong2021/100`: KRR prediction results for Kong2021 in HCP dataset with 100 resolution (3 splits).
+ `<output_dir>/KRR/HCP/Kong2021/200`: KRR prediction results for Kong2021 in HCP dataset with 200 resolution (3 splits).
+ `<output_dir>/KRR/ABCD/PrincipalGrad/1`: KRR prediction results for PrincipalGrad in ABCD dataset with 1 gradient.
+ `<output_dir>/LRR/HCP/Schaefer2018/100`: LRR fracridge prediction results for Schaefer2018 in HCP dataset with 100 resolution (3 splits).
+ `<output_dir>/LRR/HCP/Schaefer2018/200`: LRR fracridge prediction results for Schaefer2018 in HCP dataset with 200 resolution (3 splits).
+ `<output_dir>/LRR/ABCD/sICA/50`: LRR fracridge prediction results for sICA in ABCD dataset with 50 resolution.
+ `<output_dir>/LRR/ABCD/sICA/100`: LRR fracridge prediction results for sICA in ABCD dataset with 100 resolution.
+ `<output_dir>/KRR/Kong2021_opt_res/1`: KRR prediction results for Kong2021 in HCP dataset with optimized resolution (for the first split).
+ `<output_dir>/LRR/sICA_opt_res/1`: LRR fracridge prediction results for sICA in ABCD dataset with optimized resolution.

----

Check Results
====

The reference results were saved in the CBIG server: 
`$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/ref_results`

----

Bugs and Questions
====
Please contact Ru(by) Kong at roo.cone@gmail.com.

