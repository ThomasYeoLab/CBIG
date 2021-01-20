## Reference
Siyi Tang*, Nanbo Sun*, Dorothea L. Floris, Xiuming Zhang, Adriana Di Martino, B.T. Thomas Yeo. [Reconciling Dimensional and Categorical Models of Autism Heterogeneity: a Brain Connectomics & Behavioral Study](https://doi.org/10.1016/j.biopsych.2019.11.009). Biological psychiatry, in press.

----

## Background
Autism spectrum disorder (ASD) is characterized by abnormal brain connectivity, social impairment and repetitive behaviors. The range and severity of symptoms, and the pattern of abnormal brain connectivity can vary widely across autistic individuals. This variability has hindered the development of effective biomarkers. We applied a new mathematical model to functional brain imaging data of a large international cohort of autistic individuals, revealing three distinct patterns of abnormal brain connectivity. The different connectivity patterns explained variation in social impairment, repetitive behavior, executive dysfunction (inability to plan and execute plans), affective problems and externalizing problems across ASD individuals.

----

## Data Release

Factor compositions of the ASD participants included in this study are released as a spreadsheet in `data_release/factorCompositions_ASD.csv`. Specifically, the spreadsheet includes factor compositions (i.e., Pr(Factor|Participant)) of the 306 ASD participants in ABIDE-II+GENDAAR datasets as well as 166 ASD participants in ABIDE-I dataset.

The factor-specific hypo/hyper resting-state functional connectivity (RSFC) patterns (i.e., E(RSFC patterns|Factor)) and the estimated model parameters in our study are released in the folder `data_release/files_polarLDA`.

----

## Code Release
The code utilized in this study includes the following three steps. Please see `README.md` in each folder for details on each step.

### Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: https://github.com/ThomasYeoLab/Standalone_Tang2020_ASDFactors

### Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

To download the version of the code that was last tested, you can either

* visit this link:  https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.18.1-Update_stable_project_unit_test

or

* run the following command, if you have Git installed
```
git checkout -b Tang2020_ASDFactors v0.18.1-Update_stable_project_unit_test
```

### Usage
#### Step 1: Converting RSFC data to "documents"
* Folder `step1_FC2doc` contains the functions to convert RSFC data to "documents" for factor estimate or factor composition inference in next step.

#### Step 2: Estimating latent factors, visualizing factors & inferring factor compositions of new participants
* Folder `step2_FC2doc` contains the functions to:
  1) Estimate latent factors using the Bayesian model proposed in our study;
  2) Visualize the latent factors;
  3) Infer factor compositions of new ASD participants.
* The source code of the Bayesian model proposed in our study is released in `$CBIG_CODE_DIR/external_packages/polarlda-c-dist/code`.

#### Step 3: Post analyses
* Folder `step3_analyses` contains post analyses in our study. This folder is mainly for CBIG lab internal use.

#### Examples
In the folder `examples`, we provide toy example code on how to convert RSFC data to "documents", estimate & visualize latent factors, as well as infer factor compositions of unseen ASD participants.

----

## Updates
* Release v0.15.0 (07/10/2019): Initial release of Tang2020_ASDFactors
* Release v0.15.2 (13/10/2019): Update README and figure utils
* Release v0.17.0 (19/02/2020): Avoid using absolute paths. Add new environment variables to avoid possible problems caused by hard-coded absolute paths.
* Release v0.18.1 (20/01/2021): Update unit test to accommodate to the new HPC.
----

## Bugs and Questions
Please contact Siyi Tang at tangsy935@gmail.com, Nanbo Sun at sun464879934@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
