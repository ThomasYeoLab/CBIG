This folder contains the code to run variance component model on each dataset separately. 

Note that in our paper, a jackknife procedure was utilized to quantify the uncertainty of the variance estimates. For each trait (i.e. the target variable `y`, including behavioral measures, age, or sex in our paper), half the participants were randomly removed. The variance component model was fitted for the remaining participants. This jackknife procedure was repeated 1000 times, resulting in 1000 jackknife estimates of the explained trait variance.

# Code

## C1. GSP scripts

- `scripts/CBIG_LiGSR_LME_workflowGSP.sh`: 

  the top-level wrapper script for the **GSP** dataset, calling `scripts/CBIG_LiGSR_LME_workflowGSP.m`. Type `CBIG_LiGSR_LME_workflowGSP.sh` in command line for the details of usage.

- `scripts/CBIG_LiGSR_LME_workflowGSP.m`: 

  the matlab script that integrates computing functional similarity matrix (FSM), generating jackknife samples, and performing variance component model on each jackknife sample for the **GSP** dataset.

- `scripts/CBIG_LiGSR_explained_variance_GSP.m`: 

  given a subject list and trait list, this script performs variance component model specifically in the **GSP** dataset.

## C2. HCP scripts

- `scripts/CBIG_LiGSR_LME_workflowHCP.sh`: 

  the top-level wrapper script for the **HCP** dataset, calling `scripts/CBIG_LiGSR_LME_workflowHCP.m`. Type `CBIG_LiGSR_LME_workflowHCP.sh` in command line for the details of usage.

- `scripts/CBIG_LiGSR_LME_workflowHCP.m`: 

  the matlab script that integrates computing FSM, generating jackknife samples, and performing variance component model on each jackknife sample for the **HCP** dataset.

- `scripts/CBIG_LiGSR_explained_variance_HCP.m`: 

  given a subject list and trait list, this script performs variance component model specifically in the **HCP** dataset.

## C3. Scripts shared across datasets

- `scripts/utilities/CBIG_LiGSR_compute_FSM_from_FC.m`: 

  given a resting-state functional connectivity (RSFC) file, this script computes the FSM matrix based on the lower-triangular, off-diagonal entries.

- `scripts/utilities/CBIG_LiGSR_NchooseD_families.m`: 

  given a subject list, this script generates the list of subjects to be removed for each jackknife sample, considering family structure.

- `scripts/utilities/CBIG_LiGSR_LME_cmp2pipe_allstats.m`: 

  this script calculates the statistics mentioned in Li et al. (under review) related to variance component model.

- `scripts/utilities/CBIG_LiGSR_del_d_jack_cmp2pipelines.m`: 

  it computes 
  1. the percentage improvement of the explained variance (baseline+GSR pipeline versus baseline pipeline)
  2. the jackknife mean and variance of improved explained variance in (i). 
  
  This script is called by `CBIG_LiGSR_LME_cmp2pipe_allstats.m`.

- `scripts/utilities/CBIG_LiGSR_PosNeg_jackIQR_cmp2pipe.m`: 

  it computes (1) the number of traits with entire interquartile range (IQR) of explained variance **difference** above 0 (or below 0); (2) the number of traits with the median of explained variance **difference** above 0 (or below 0).
  
# Usage

## U1. Run variance component analysis using the GSP dataset

### U1.1 Preparations
1. The RSFC files of both preprocessing pipelines (baseline, baseline+GSR) need to be precomputed. For each pipeline, the RSFC file is a `.mat` file containing a 3-D matrix with the dimensions of RxRxN, where R represents the number of ROIs and N is the number of subjects selected by the user.
2. The CSV file containing the demographic and behavioral information should be downloaded from the GSP website.
3. The list of subject IDs the user would like to include needs to be prepared in the form of a text file. In this list, each line corresponds to one subject ID.
4. The trait names need to be written in a text file, where each line corresponds to one trait name. Note that the trait names are obtained from the headers of the CSV file the user downloaded from the GSP website. For example, `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/GSP_lists/intelligence_score.txt`.
5. The covariate names need to be written in a separate text file, where each line corresponds to one covariate name. Note that the covariate names are obtained from the headers of the CSV file the users downloaded from the GSP website. If it is necessary to consider framewise displacement (FD) as a covariate, you have to add a line "FD" in this text file. If it is necessary to consider DVARS as a covariate, you have to add a line "DVARS" in this text file. An example of this text file is `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/scripts/GSP_lists/covariates.txt`.
6. If the covariates include FD, the mean FD of individual subjects need to be written as a column in a separate text file (i.e. each line corresponds to the mean FD of a subject). Same for the mean DVARS of each subject.

### U1.2 Run the GSP workflow
Use the following command for each preprocessing pipeline

```
$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowGSP.sh \
  -RSFC_file <RSFC_file> -trait_list <trait_list> -covariate_list <covariate_list> \
  -FD_file <FD_file> -DVARS_file <DVARS_file> -subject_list <subject_list> -outdir <outdir> -d <d> \
  -num_samples <num_samples> -rmsub_prefix <rmsub_prefix> -data_csv <data_csv> -ystem <ystem>
```

- `RSFC_file`: Absolute path to the RSFC (`.mat`) file generated using the pipeline (baseline or baseline+GSR).
- `triat_list`: Absolute path to the text file containing the list of trait names.
- `covariate_list`: Absolute path to the text file containing the list of covariate names.
- `FD_file`: Absolute path to the text file containing mean FD of each subject, if FD needs to be included as a covariate. Pass in "NONE" if the user does not want to include FD as a covariate.
- `DVARS_file`: Absolute path to the text file containing mean DVARS of each subject, if DVARS needs to be included as a covariate. Pass in "NONE" if the user does not want to include DVARS as a covariate.
- `subject_list`: Absolute path to the text file containing the list of subject IDs.
- `outdir`: Absolute path of output directory. It is recommended for the user to pass in different output directories for different pipelines.
- `d`: The number of subjects to be removed for each jackknife sample, e.g. 431.
- `num_samples`: The total number of jackknife samples, e.g. 1000.
- `rmsub_prefix`: The prefix to the filename of the subject list to be removed for jackknife samples. The list of removed subject IDs of each jackknife sample will be saved under the output file name `${outdir}/jackknife_lists/${rmsub_prefix}_choose${d}_set${i}.txt`, where `i` ranges from 1 to `num_samples`.
- `data_csv`: Absolute path to the CSV file the user downloaded from the GSP website.
- `y_stem`: The trait values will be read from `data_csv` and saved in a .mat file: `${outdir}/y_${ystem}.mat`. For example, if `trait_list` contains 23 behavioral names, you can set `ystem = 23behaviors`.

For each preprocessing pipeline, this command will create the following outputs inside `outdir`:

- A folder called `fullset`, which contains M `.mat` files (M is the total number of traits). Each `.mat` file stores the explained variance for a trait estimated based on the full set of subjects.
- `num_samples` folders called `del${d}_set${i}`, where `i` ranges from 1 to `num_samples`. Each folder contains M `.mat` files. Each `.mat` file stores the explained variance for a trait estimated on the `i`-th jackknife sample.
- A folder called `jackknife_lists`, which contains `num_sample` text files. Each text file is a list of subjects that were removed from the full set for the corresponding jackknife sample.
- A folder called `logs`, which contains the echoed information from the `CBIG_LiGSR_LME_workflowGSP.sh` command.
- `covariates_${covariate_names}.mat`, where `covariate_names` is a string of the concatenated covariate names, e.g. `Age_Bin_DVARS_FD_Sex`. This `.mat` file contains the covariate matrix that were plugged into the variance component model.
- `FSM.mat`, which is the functional similarity matrix among all subjects.
- `y_${y_stem}.mat`, which is the values of all traits extracted from the GSP CSV file.

### U1.3 Compare explained variances between preprocessing pipelines
Run the following command in matlab terminal

```
cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities'));
CBIG_LiGSR_LME_cmp2pipe_allstats( trait_list, pipeline1_dir, pipeline2_dir, num_families, num_samples, d, outdir );
```

- `trait_list`: A string, the absolute path to a text file containing the list of trait names.
- `pipeline1_dir`: A string, the `outdir` variable the user passed to the `CBIG_LiGSR_LME_workflowGSP.sh` command for the first preprocessing pipeline (i.e. baseline+GSR).
- `pipeline2_dir`: A string, the `outdir` variable the user passed to the `CBIG_LiGSR_LME_workflowGSP.sh` command for the second preprocessing pipeline (i.e. baseline).
- `num_families`: A scalar, the number of families in the full set of subjects. Since the subjects in the GSP data are unrelated, the number of families is the same as total number of subjects.
- `num_samples`: A scalar, the total number of jackknife samples (a scalar), e.g. 1000.
- `d`: A scalar, the number of subjects to be removed for each jackknife sample (a scalar), e.g. 431.
- `outdir`: A string, the full-path folder name to save the output `allstats_cmp2pipelines.mat` file. This `.mat` file contains the following variables
  - `mean_m2`: the mean explained variance averaged across all traits and all jackknife samples for each preprocessing pipeline;
  - `m_jack`: the jackknife mean of the difference in the explained variance between the two preprocessing pipelines;
  - `v_jack`: the jackknife variance of the difference in the explained variance between the two preprocessing pipelines;
  - `perc_improv`: the percentage improvement of the mean explained variance using pipeline 1 compared to pipeline 2;
  - `IQR_pos` and `IQR_neg`: the number of traits with entire interquartile range (IQR) of explained variance **difference** above 0 (or below 0);
  - `med_pos` and `IQR_neg`: the number of traits with the median of explained variance **difference** above 0 (or below 0)

## U2. Run variance component analysis using the HCP dataset

### U2.1 Preparations
1. The RSFC files of both preprocessing pipelines (baseline, baseline+GSR) need to be precomputed. For each pipeline, the RSFC file is a `.mat` file containing a 3-D matrix with the dimensions of RxRxN, where R represents the number of ROIs and N is the number of subjects selected by the user.
2. The CSV files (both the restricted and unrestricted ones) containing the demographic and behavioral information should be downloaded from the HCP website.
3. The list of subject IDs the user would like to include needs to be prepared in the form of a text file. In this list, each line corresponds to one subject ID. It is suggested to select only unrelated subjects because dealing with both RSFC and family structure within the variance component model is tricky.
4. The trait names need to be written in a text file, where each line corresponds to one trait name. Note that the trait names are obtained from the headers of the CSV files the user downloaded from the HCP website. For example, `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/replication/scripts/HCP_lists/Cognitive_unrestricted.txt`.
5. The covariate names need to be written in a separate text file, where each line corresponds to one covariate name. Note that the covariate names are obtained from the headers of the CSV files the users downloaded from the HCP website. If it is necessary to consider framewise displacement (FD) as a covariate, you have to add a line "FD" in this text file. If it is necessary to consider DVARS as a covariate, you have to add a line "DVARS" in this text file. An example of this text file is `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/replication/scripts/HCP_lists/covariates_58behaviors.txt`.
6. If the covariates include FD, the mean FD of individual subjects need to be written as a column in a separate text file (i.e. each line corresponds to the mean FD of a subject). Same for the mean DVARS of each subject.

### U2.2 Run the HCP workflow
Use the following command for each preprocessing pipeline

```
$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/VarianceComponentModel/scripts/CBIG_LiGSR_LME_workflowHCP.sh \
  -RSFC_file <RSFC_file> -trait_list <trait_list> -covariate_list <covariate_list> -FD_file <FD_file> \
  -DVARS_file <DVARS_file> -subject_list <subject_list> -outdir <outdir> -d <d> -num_samples <num_samples> \
  -rmsub_prefix <rmsub_prefix> -restricted_csv <restricted_csv> -unrestricted_csv <unrestricted_csv> \
  -unrelated <unrelated> -ystem <ystem>
```

- `RSFC_file`: Absolute path to the RSFC file (`.mat`) generated using the pipeline (baseline or baseline+GSR).
- `trait_list`: Absolute path to the text file containing the list of trait names.
- `covariate_list`: Absolute path to the text file containing the list of covariate names.
- `FD_file`: Absolute path to the text file containing mean FD of each subject, if FD needs to be included as a covariate. Pass in "NONE" if the user does not want to include FD as a covariate.
- `DVARS_file`: Absolute path to the text file containing mean DVARS of each subject, if DVARS needs to be included as a covariate. Pass in "NONE" if the user does not want to include DVARS as a covariate.
- `subject_list`: Absolute path to the text file containing the list of subject IDs.
- `outdir`: Absolute path of output directory. It is recommended for the user to pass in different output directories for different pipelines.
- `d`: The number of subjects to be removed for each jackknife sample, e.g. 209.
- `rmsub_prefix`: The prefix to the filename of the subject list to be removed for jackknife samples. The list of removed subject IDs of each jackknife sample will be saved under the output file name `${outdir}/jackknife_lists/${rmsub_prefix}_choose${d}_set${i}.txt`, where `i` ranges from 1 to `num_samples`.
- `restricted_csv`: Absolute path to the restricted CSV file the user downloaded from the HCP website.
- `unrestricted_csv`: Absolute path to the unrestricted CSV file the users downloaded from the HCP website.
- `unrelated`: 0 or 1. 1 means all subjects in the `subject_list` are unrelated. 0 means there is a family structure within the subjects in `subject_list`, and the family information will be read from `restricted_csv`.
- `y_stem`: The trait values will be read from `restricted_csv` and `unrestricted_csv`, and saved in a .mat file: `${outdir}/y_${ystem}.mat`. For example, if '"'trait_list'"' contains 13 cognitive behavioral names, you can set `ystem = 13cognitive`.

For each preprocessing pipeline, this command will create the following outputs inside `outdir`:

- A folder called `fullset`, which contains M `.mat` files (M is the total number of traits). Each `.mat` file stores the explained variance for a trait estimated based on the full set of subjects.
- `num_samples` folders called `del${d}_set${i}`, where `i` ranges from 1 to `num_samples`. Each folder contains M `.mat` files. Each `.mat` file stores the explained variance for a trait estimated on the `i`-th jackknife sample.
- A folder called `jackknife_lists`, which contains `num_sample` text files. Each text file is a list of subjects that were removed from the full set for the corresponding jackknife sample.
- A folder called `logs`, which contains the echoed information from the `CBIG_LiGSR_LME_workflowHCP.sh` command.
- `covariates_${covariate_names}.mat`, where `covariate_names` is a string of the concatenated covariate names, e.g. `Age_in_Yrs_DVARS_FD_Gender`. This `.mat` file contains the covariate matrix that were plugged into the variance component model.
- `FSM.mat`, which is the functional similarity matrix among all subjects.
- `y_${y_stem}.mat`, which is the values of all traits extracted from the HCP CSV files.

### U2.3 Compare explained variances between preprocessing pipelines
Run the following command in matlab terminal

```
cd(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities'));
CBIG_LiGSR_LME_cmp2pipe_allstats( trait_list, pipeline1_dir, pipeline2_dir, num_families, num_samples, d, outdir );
```

- `trait_list`: A string, the absolute path to a text file containing the list of trait names (a string).
- `pipeline1_dir`: A string, the `outdir` variable the user passed to the `CBIG_LiGSR_LME_workflowHCP.sh` command for the first preprocessing pipeline (i.e. baseline+GSR).
- `pipeline2_dir`: A string, the `outdir` variable the user passed to the `CBIG_LiGSR_LME_workflowHCP.sh` command for the second preprocessing pipeline (i.e. baseline).
- `num_families`: A scalar, the number of families in the full set of subjects. 
- `num_samples`: A scalar, the total number of jackknife samples (a scalar), e.g. 1000.
- `d`: A scalar, the number of subjects to be removed for each jackknife sample (a scalar), e.g. 209.
- `outdir`: A string, the full-path folder name to save the output `allstats_cmp2pipelines.mat` file. This `.mat` file contains the following variables
  - `mean_m2`: the mean explained variance averaged across all traits and all jackknife samples for each preprocessing pipeline;
  - `m_jack`: the jackknife mean of the difference in the explained variance between the two preprocessing pipelines;
  - `v_jack`: the jackknife variance of the difference in the explained variance between the two preprocessing pipelines;
  - `perc_improv`: the percentage improvement of the mean explained variance using pipeline 1 compared to pipeline 2;
  - `IQR_pos` and `IQR_neg`: the number of traits with entire interquartile range (IQR) of explained variance **difference** above 0 (or below 0);
  - `med_pos` and `IQR_neg`: the number of traits with the median of explained variance **difference** above 0 (or below 0)
  
  
  
