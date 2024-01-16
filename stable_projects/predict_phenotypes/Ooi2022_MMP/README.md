# Comparison of individualized behavioral predictions across anatomical, diffusion and functional connectivity MRI (MMP)
# REFERENCE
* Ooi, L.Q.R., Chen, J., Zhang, S., Tam, A., Li, J., Dhamala, E., Zhou, J.H., Holmes, A.J., Yeo, B.T. [Comparison of individualized behavioral predictions across anatomical, diffusion and functional connectivity MRI](https://doi.org/10.1016/j.neuroimage.2022.119636). Neuroimage, 2022. 

# BACKGROUND
We extracted features derived from anatomical T1, diffusion and functional imaging, and compared their ability in predicting behaviour at an individual level using the HCP and ABCD using 3 regression models. We found that functional MRI gave the best prediction results, regardless of type of regression or behaviour. The combination of these modalities through stacking improved predictions of cognition, but not other aspects of behaviour. Additionally, combining features of functional MRI performs as well as combining all features from all modalities. 

# CODE RELEASE
## Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: [https://github.com/ThomasYeoLab/Standalone_Ooi2022_MMP](https://github.com/ThomasYeoLab/Standalone_Ooi2022_MMP)

## Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

To download the version of the code that was last tested, you can either visit this link: 

[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.29.6-He2022_Ooi2022-updates](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.29.6-He2022_Ooi2022-updates)

run the following command, if you have Git installed

```
git checkout -b He2022_Ooi2022-updates v0.29.6-He2022_Ooi2022-updates
```

# USAGE
## Setup
1. Make sure you have installed: Matlab 2018b
2. Follow `$CBIG_CODE_DIR/setup/README.md` to setup the config file. Instead of using `CBIG/setup/CBIG_sample_config.sh`, you need to use `$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Ooi2022_MMP/replication/config/CBIG_MMP_tested_config.sh`.
3. The data is stored in the following hyperlinks for the [HCP](https://github.com/ThomasYeoLab/Ooi2022_MMP_HCP) and [ABCD](https://dx.doi.org/10.15154/1523482). The HCP repository is publicly available. However, please apply for restricted access to the HCP to gain access to information such as age and the Family ID. The ABCD data is available on the NDA website through in study [1368](https://dx.doi.org/10.15154/1523482).

## Replication
* `replication` this folder contains scripts to replicate all the analysis in this paper. This requires the user to have access to the ABCD and HCP repositories.

## Examples and unit tests
* `examples` this folder provides a toy example for running all regression codes.
* `unit_tests` this folder runs codes in `examples` and checks the reference output.

## Regression models
* `regression`: this folder contains the workflows for KRR, LRR, Elasticnet, multiKRR and stacking for both the HCP and ABCD. Additionally, it contains the scripts to perform permutation tests for KRR, multiKRR and stacking models to verify if they performed above chance.
## Generate figures 
* `generate_figures` this folder contains the scripts to generate the figures in our paper given you have all the results


# UPDATES
- Release v0.29.6 (16/1/2024): Addition of code for feature importance plotting
- Release v0.24.0 (23/9/2022): Initial release of Ooi2022_MMP project

# BUGS and QUESTIONS
Please contact Leon Ooi at leonooiqr@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
