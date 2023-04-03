# Examples of comparison between gradients and parcellations for RSFC behavioral prediction

This example shows how to generate prediction results for a given parcellation with two different resolutions 100 and 200. The example also shows how to optimize the resolution for this parcellation. In this example, we will perform behavioral prediction for 30 example particiants. There are 3 target variables to be predicted. The prediction will be performed using KRR and LRR with 3-fold cross-validation. The prediction results will be saved in the `output` folder.

----

References
==========
+ Kong et al., [Comparison Between Gradients and Parcellations for Functional Connectivity Prediction of Behavior](https://doi.org/10.1016/j.neuroimage.2023.120044). NeuroImage. 2023

----

Data
====
We prepared the following data for 30 participants:
+ feature matrices (the lower triangular part of FC matrcies with 100 and 200 ROIs) for each participant.
+ target variables (3 behavioral measures) for each participant.
+ 3-fold cross-validation setup for 30 participants.
+ covariate matrices for each participant (4 regressors).

----

Run
====

### Perform prediction and optimize resolution
----

The following script will generate prediction results for 100 and 200 resolutions for KRR and LRR. After that, the script will optimize the resolution for this parcellation. Check the script for more details about how we setup the parameters.


```
In the terminal, specify the output directory and call the script:

```
$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/examples/CBIG_GradPar_example_KRR_wrapper.m <output_dir>

$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/examples/CBIG_GradPar_example_LRR_frac_wrapper.m <output_dir>
```

### Check example results
----

We provide a script which checks user's example results with our reference results:

```
$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Kong2023_GradPar/examples/CBIG_GradPar_check_example_results.m <output_dir>
```

----

Bugs and Questions
====
Please contact Ru(by) Kong at roo.cone@gmail.com.

