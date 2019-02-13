References
==========
+ Li J, Kong R, Liegeois R, Orban C, Tan Y, Sun N, Holmes AJ, Sabuncu MR, Ge T, Yeo BTT. **Global Signal Regression Strengthens Association between Resting-State Functional Connectivity and Behavior**, under review.

----

Background
==========
Global signal regression is a controversial preprocessing step for resting-state functional MRI. GSR effectively removes global motion and respiratory artifacts. However, GSR might distort resting-state functional connectivity (RSFC) by introducing negative correlations and removing neural information. The vast majority of studies have focused on the effectiveness of GSR in removing imaging artifacts, and its potential biases. Given growing interest in functional connectivity fingerprinting, here we explored the utilitarian question of whether GSR strengthens or weakens the association between RSFC and various behavioral measures. Two large-scale datasets were involved: the Brain Genomics Superstruct Project (GSP; N=862; 23 behavioral measures) and the Human Connectome Project (HCP; N=953; 58 behavioral measures). Using variance component model and kernel ridge regression, we showed that GSR enhanced the associations between whole-brain RSFC and behaviors as well as the prediction accuracies of many behavioral measures, at least in the case for young healthy subjects. Since GSR was more effective at removing motion-related and respiratory-related artifacts, these improvements were unlikely to be the result of motion-related or respiratory-related artifacts.

----

Code Release
===========
- `VarianceComponentModel/scripts` folder contains codes to run variance component model on each dataset separately. Specifically, `CBIG_LiGSR_LME_workflowGSP.sh` and `CBIG_LiGSR_LME_workflowHCP.sh` are the top-level wrappers calling other scripts.

- `KernelRidgeRegression` folder contains codes to run kernel ridge regression to predict behavioral measures for each dataset. `KernelRidgeRegression/GSP/scripts/CBIG_LiGSR_KRR_workflowGSP.sh` and `KernelRidgeRegression/HCP/scripts/CBIG_LiGSR_KRR_workflowHCP.sh` are the top-level wrappers.

- `unit_tests` folder contains the codes for two types of unit tests: (1) replication of results shown in Li et al. (under review); (2) results of 1 or 2 behavioral measures related to intelligence for each dataset. For the procedures to run the replication tests, refer to `unit_tests/replication/README.md`. For the procedures to run the tests on only the 1 or 2 intelligence measures, refer to `unit_tests/intelligence_score/README.md`.

- `config` folder contains the latest configuration files and matlab startup.m that were successfully tested.

Note that this project uses generic functions from other folders, which may be updated over time. To download the version of the code that was last tested, you can either

- visit this link: https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.9.0-Li2019_GSR

or

- run the following command, if you have Git installed
```
git checkout -b Li2019_GSR v0.9.0-Li2019_GSR
```

----

Updates
===========

- Release v0.9.0 (13/02/2019): Initial release of Li2019_GSR project.

Bugs and Questions
==========
Please contact Jingwei Li at jingweili.sjtu.nus@gmail.com, Ru Kong at roo.cone@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
