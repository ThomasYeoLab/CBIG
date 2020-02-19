This README includes the instruction on how to run the unit tests for Li2019_GSR project for 1 or 2 behavioral measures related to intelligence for each dataset. This is a faster version of the unit tests in `../replication` folder.

Note that due to the restriction of the datasets as well as the size of input files, only users that have access to our servers can run the unit tests.

This unit test needs a pre-defined envrionmental variance `$CBIG_TESTDATA_DIR`. It has been specified in the config file: `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/config/CBIG_LiGSR_tested_config.sh`.

## Variance component model

This section shows how to estimate the explained variances of the 1 or 2 measures related to intelligence using variance component model.

### The Brain Genomics Superstruct Project (GSP)
#### Data
This unit test will estimate the RSFC-explained variance of two measures: Shipley Vocabulary and WAIS - Matrix Reasoning. 

It assumes that the user has obtained the CSV file with all behavioral and demographic information from the GSP dataset. An example file (shuffled and noise added):

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/VarianceComponentModel/GSP/lists/GSP_phenotypes_shuffled.csv`

It is also assumed that the resting-state functional connectivity matrix was pre-computed:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/VarianceComponentModel/GSP/RSFC_862_Fisher_Baseline.mat` using the baseline preprocessing;

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/VarianceComponentModel/GSP/RSFC_862_Fisher_GSR.mat` using the preprocessing pipeline with GSR.

#### Code
- Step 1. Perform variance component model for each preprocessing pipeline
  
  If you have the access to our server (because the below stated function submits jobs to circ-spool, i.e. our HPC cluster), run the following command in the terminal:
  
  `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/CBIG_LiGSR_LME_unittest_intelligence_score_GSP.sh   <output_dir>`
  
  where `<output_dir>` is the output directory that the user can specify.

- Step 2. Integrate and compare the results between the two pipelines
  
  Execute this step after all jobs in Step 1 are finished.
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LME_unittest_intelligence_score_cmp2pipe_GSP( output_dir )`
  
  where `output_dir` is the output directory the user specified in Step 1.

- Step 3. Compare your results with the groud truth results
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LME_unittest_intelligence_score_cmp_w_reference_GSP( all_stats_mat )`
  
  where `all_stats_mat` is the .mat file generated in Step 2. It should be stored as `[output_dir '/compare_2pipe/allstats_cmp2pipelines.mat']`. A sucessful message `Your results replicated the reference results.` will be printed out if you generated the correct explained variances.
  

### The Human Connectome Project (HCP)
#### Data
This unit test will estimate the explained variance of fluid intelligence measure in the HCP dataset (PMAT24_A_CR).

It assumes that the user has two CSV files containing the phenotypes and are stored in the following locations:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/VarianceComponentModel/HCP/lists/fake_family.csv`

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/VarianceComponentModel/HCP/lists/fake_behavior.csv`

It is also assumed that the resting-state functional connectivity matrix was pre-computed:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/replication/VarianceComponentModel/HCP/RSFC_953_unrelated_419_Fisher_Baseline.mat` using the baseline preprocessing;

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/replication/VarianceComponentModel/HCP/RSFC_953_unrelated_419_Fisher_GSR.mat` using the preprocessing pipeline with GSR.

#### Code
- Step 1. Perform variance component model for each preprocessing pipeline
  
  If you have the access to our server (because the below stated function submits jobs to circ-spool, i.e. our HPC cluster), run the following command in the terminal:
  
  `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/CBIG_LiGSR_LME_unittest_PMAT_HCP.sh   <output_dir>`
  
  where `<output_dir>` is the output directory that the user can specify.

- Step 2. Integrate and compare the results between the two pipelines
  
  Execute this step after all jobs in Step 1 are finished.
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LME_unittest_PMAT_cmp2pipe_HCP( output_dir )`
  
  where `output_dir` is the output directory the user specified in Step 1.
  
- Step 3. Compare your results with the groud truth results
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LME_unittest_PMAT_cmp_w_reference_HCP( all_stats_mat )`
  
  where `all_stats_mat` is the .mat file generated in Step 2. It should be stored as `[output_dir '/compare_2pipe/allstats_cmp2pipelines.mat']`. A sucessful message `Your results replicated the reference results.` will be printed out if you generated the correct explained variances.
  

## Kernel ridge regression
The section below will state the steps required to predict the 1 or 2 measures related to intelligence using kernel ridge regression in two different datasets.

Please pass in different `output_dir` (see below) for different datasets to avoid overwrites.

### The Brain Genomics Superstruct Project (GSP)
#### Data
This unit test will predict two measures: Shipley Vocabulary and WAIS - Matrix Reasoning. 

It assumes that the user has obtained the CSV file with all behavioral and demographic information from the GSP dataset. An example file (shuffled and noise added):

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/GSP/lists/GSP_phenotypes_shuffled.csv`

It is also assumed that the resting-state functional connectivity matrix was pre-computed:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/GSP/RSFC_862_Fisher_Baseline.mat` using the baseline preprocessing;

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/GSP/RSFC_862_Fisher_GSR.mat` using the preprocessing pipeline with GSR.

#### Code
- Step 1. Perform kernel ridge regression using the resting-state functional connectivity matrix generated using each of the preprocessing pipeline
  
  If you have the access to our server (because the below stated function submits jobs to circ-spool, i.e. our HPC cluster), run the following command in the terminal:
  
  `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/CBIG_LiGSR_KRR_unittest_intelligence_score_GSP.sh   <output_dir>`
  
  where `<output_dir>` is the output directory that the user can specify.

- Step 2. Integrate and compare the results between the two pipelines
  
  Execute this step after all jobs in Step 1 are finished.
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_KRR_unittest_intelligence_score_cmp2pipe_GSP( output_dir )`
  
  where `output_dir` is the output directory the user specified in Step 1.
  
- Step 3. Compare your results with the groud truth results
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_KRR_unittest_intelligence_score_cmp_w_reference_GSP( final_result )`
  
  where `final_result` is the .mat file generated in Step 2. It should be stored as `[output_dir '/compare_2pipe/final_result.mat']`. A sucessful message `Your results replicated the reference results.` will be printed out if you generated the correct kernel regression accuracies.


### The Human Connectome Project (HCP)
#### Data
This unit test will predict the fluid intelligence measure in the HCP dataset (PMAT24_A_CR).

It assumes that the user has two CSV files containing the phenotypes, e.g.

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/HCP/lists/fake_family.csv`

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/HCP/lists/fake_behavior.csv`

It is also assumed that the resting-state functional connectivity matrix was pre-computed:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/HCP/cort+subcort_new_S1200_953_Fisher_Baseline.mat` using the baseline preprocessing;

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/KernelRidgeRegression/HCP/cort+subcort_new_S1200_953_Fisher_GSR.mat` using the preprocessing pipeline with GSR.


#### Code
- Step 1. Perform kernel ridge regression using the resting-state functional connectivity matrix generated using each of the preprocessing pipeline
  
  If you have the access to our server (because the below stated function submits jobs to circ-spool, i.e. our HPC cluster), run the following command in the terminal:
  
  `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/CBIG_LiGSR_KRR_unittest_PMAT_HCP.sh   <output_dir>`
  
  where `<output_dir>` is the output directory that the user can specify.

- Step 2. Integrate and compare the results between the two pipelines
  
  Execute this step after all jobs in Step 1 are finished.
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_KRR_unittest_PMAT_cmp2pipe_HCP( output_dir )`
  
  where `output_dir` is the output directory the user specified in Step 1.
  
- Step 3. Compare your results with the groud truth results
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_KRR_unittest_PMAT_cmp_w_reference_HCP( final_result )`
  
  where `final_result` is the .mat file generated in Step 2. It should be stored as `[output_dir '/compare_2pipe/final_result.mat']`. A sucessful message `Your results replicated the reference results.` will be printed out if you generated the correct kernel regression accuracies.


## Linear ridge regression
The section below will state the steps required to predict the 1 or 2 measures related to intelligence using linear ridge regression in two different datasets.

Please pass in different `output_dir` (see below) for different datasets to avoid overwrites.

### The Brain Genomics Superstruct Project (GSP)
#### Data
This unit test will predict two measures: Shipley Vocabulary and WAIS - Matrix Reasoning. 

It assumes that the user has obtained the CSV file with all behavioral and demographic information from the GSP dataset. An example file (shuffled and noise added):

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/GSP/lists/GSP_phenotypes_shuffled.csv`

It is also assumed that the resting-state functional connectivity matrix was pre-computed:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/GSP/RSFC_862_Fisher_Baseline.mat` using the baseline preprocessing;

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/GSP/RSFC_862_Fisher_GSR.mat` using the preprocessing pipeline with GSR.

#### Code
- Step 1. Perform linear ridge regression using the resting-state functional connectivity matrix generated using each of the preprocessing pipeline
  
  If you have the access to our server (because the below stated function submits jobs to circ-spool, i.e. our HPC cluster), run the following command in the terminal:
  
  `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/CBIG_LiGSR_LRR_unittest_intelligence_score_GSP.sh   <output_dir>`
  
  where `<output_dir>` is the output directory that the user can specify.

- Step 2. Integrate and compare the results between the two pipelines
  
  Execute this step after all jobs in Step 1 are finished.
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LRR_unittest_intelligence_score_cmp2pipe_GSP( output_dir )`
  
  where `output_dir` is the output directory the user specified in Step 1.
  
- Step 3. Compare your results with the groud truth results
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LRR_unittest_intelligence_score_cmp_w_reference_GSP( final_result )`
  
  where `final_result` is the .mat file generated in Step 2. It should be stored as `[output_dir '/compare_2pipe/final_result.mat']`. A sucessful message `Your results replicated the reference results.` will be printed out if you generated the correct linear regression accuracies.


### The Human Connectome Project (HCP)
#### Data
This unit test will predict the fluid intelligence measure in the HCP dataset (PMAT24_A_CR).

It assumes that the user has two CSV files containing the phenotypes, e.g.

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/HCP/lists/fake_family.csv`

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/HCP/lists/fake_behavior.csv`

It is also assumed that the resting-state functional connectivity matrix was pre-computed:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/HCP/cort+subcort_new_S1200_953_Fisher_Baseline.mat` using the baseline preprocessing;

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/Li2019_GSR/intelligence_score/LinearRidgeRegression/HCP/cort+subcort_new_S1200_953_Fisher_GSR.mat` using the preprocessing pipeline with GSR.


#### Code
- Step 1. Perform linear ridge regression using the resting-state functional connectivity matrix generated using each of the preprocessing pipeline
  
  If you have the access to our server (because the below stated function submits jobs to circ-spool, i.e. our HPC cluster), run the following command in the terminal:
  
  `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/CBIG_LiGSR_LRR_unittest_PMAT_HCP.sh   <output_dir>`
  
  where `<output_dir>` is the output directory that the user can specify.

- Step 2. Integrate and compare the results between the two pipelines
  
  Execute this step after all jobs in Step 1 are finished.
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LRR_unittest_PMAT_cmp2pipe_HCP( output_dir )`
  
  where `output_dir` is the output directory the user specified in Step 1.
  
- Step 3. Compare your results with the groud truth results
  
  In matlab command window, run the following commands
  
  `cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'))`
  
  `CBIG_LiGSR_LRR_unittest_PMAT_cmp_w_reference_HCP( final_result )`
  
  where `final_result` is the .mat file generated in Step 2. It should be stored as `[output_dir '/compare_2pipe/final_result.mat']`. A sucessful message `Your results replicated the reference results.` will be printed out if you generated the correct linear regression accuracies.





