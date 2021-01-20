# Reference

Sun N, Mormino EC, Chen J, Sabuncu MR, Yeo BTT. [**Multi-modal latent factor exploration of atrophy, cognitive and tau heterogeneity in Alzheimer's disease**](https://www.sciencedirect.com/science/article/abs/pii/S105381191930624X). NeuroImage, 201:116043, 2019

----

# Background

Alzheimer’s disease (AD) affects 10% of the elderly population. The disease remains poorly understood with no cure. The main symptom is memory loss, but other symptoms might include impaired language and executive function (ability to plan and accomplish goals; e.g., grocery shopping). The severity of behavioral symptoms and brain atrophy (gray matter loss) can vary widely across patients. This variability complicates diagnosis, treatment, and prevention. We proposed a data-driven Bayesian model (MMLDA) to identify latent factors from atrophy patterns and cognitive deficits *simultaneously*, thus exploiting the rich dimensionality within each modality. By applying our model to ADNI dataset, we found three atrophy-cognitive factors and these factors are also associated with distinct tau depositions. This model can potentially be applied to understand brain disorders with varying symptoms, including autism and schizophrenia.

----

# Data Release
Factor compositions of the ADNI subjects involved in this study are released as a spreadsheet in the `ADNIDataRelease` folder. Specifically, the spreadsheet `FactorCompositions_ADNIGO2_895bl-668m12.csv` includes 895 ADNIGO/2-enrolled subjects’ baseline factor compositions and their factor compositions 12 months after baseline (N = 668). The spreadsheet `FactorCompositions_ADNI1_777bl.csv` includes 777 ADNI1-enrolled subjects' baseline factor compositions. **Please note that the factors of ADNI1 are slightly different from ADNIGO/2, see our paper for more info**.

Folder `ADNIDataRelease/MMLDA_files/ADNIGO2/k3/` contain the probabilistic atrophy maps and cognitive deficits of the latent factors for the three-factor model. Maps for the two- and four-factor models are also released (folders `k2` and `k4`). You can also find the probabilistic atrophy maps and cognitive deficits of the three-factor model in ADNI1 here: `ADNIDataRelease/MMLDA_files/ADNI1/k3/`.

----

# Code Release
### Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: [https://github.com/ThomasYeoLab/Standalone_Sun2019_ADJointFactors](https://github.com/ThomasYeoLab/Standalone_Sun2019_ADJointFactors)

### Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

- To download the version of the code that was last tested, you can either

    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test)

    or

    - run the following command, if you have Git installed
 
    ```
    git checkout -b Sun2019_ADJointFactors v0.18.1-Update_stable_project_unit_test
    ```

----

# Usage
### Setup
1. Make sure you have installed: Matlab 2014a, Freesurfer 6.0, SPM 12, FSL 5.0.8, Ants 2.2.0
2. Install CAT12 into SPM 12 by following [CAT12 Download](http://www.neuro.uni-jena.de/cat/index.html#DOWNLOAD)
3. Follow `$CBIG_CODE_DIR/setup/README.md` to setup the config file. Instead of using `CBIG/setup/CBIG_sample_config.sh`, 
you need to use `$CBIG_CODE_DIR/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/replication/config/CBIG_MMLDA_tested_config.sh`.

### Infer factor compositions of new subjects
Even though we provide the factor compositions of ADNI2 bl/m12 and ADNI1 bl subjects, the user
may want to infer factor compositions of an unseen subject by using our estimated factors. 
The user can check `infer_new_subjects/README.md` for instructions.

### Estimate new factors
If the users want to estimate new factors or replicate our results, they can follow
4 steps.
##### step1 SPM VBM
- Folder `step1_SPM_VBM` contains the functions to do SPM VBM by using CAT12 toolbox.

##### step2 MMLDA
- Folder `step2_MMLDA` contains functions to:
    1) Convert grey matter density maps and behavioral scores to "documents"
    2) Estiamte latent factors using MMLDA 
    3) Visualize the latent factors
    4) Infer factor compositions for new participants
- The source code of the MMLDA model is released in `$CBIG_CODE_DIR/external_packages/mmlda-c-dist/code`

##### step3 PET preprocessing
- Folder `step3_PET_preprocess` contains functions to do PET preprocessing.

##### step4 Post analyses
* Folder `step4_analyses` contains post-hoc analyses in our study. This folder is mainly for CBIG lab internal use.

### Unit tests
The unit test is mainly for internal use.

### Examples
In the folder `examples`, we provide toy example code on how to estimate & visualize latent factors,
 as well as inferring factor compositions of new participants. The user can run this example to verify the set up is correct. 

----

# Updates
- Release v0.10.0 (05/04/2019): Initial release of Sun2019_ADJointFactors project
- Release v0.10.4 (07/06/2019): Update README.md
- Release v0.15.3 (16/10/2019): Update reference
- Release v0.17.0 (19/02/2020): Avoid using absolute paths. Add new environment variables to avoid possible problems caused by hard-coded absolute paths.
- Release v0.18.1 (20/01/2021): Update unit test to accommodate to the new HPC.


----

# Bugs and Questions

Please contact Nanbo Sun at sun464879934@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

