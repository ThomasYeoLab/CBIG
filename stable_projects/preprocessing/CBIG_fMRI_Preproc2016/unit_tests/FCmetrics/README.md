This document details the steps to run the preprocessing pipeline which will generate FC metrics and verify the accuracy of the code by comparing the generated FC metrices with the ground truth

Notice that all the filenames and directories below work for CBIG lab only. It is assumed that the input surface data are in `fsaverage5` space.

----

## Data

- Data set

The `CBIG_fMRI_Preproc2016` preprocessing scripts use the Brain Genomics Superstruct Project (GSP) data as a baseline data set. All data mentioned here are based on GSP data set.

GSP data contain both structual MRI (T1) and functional MRI (T2*). All subjects (N = 1570) are healthy, young subjects (age: 18-35). The preprocessed (by `recon-all`) structual MRI data and raw functional MRI are in this folder:

```
/mnt/eql/yeo1/data/GSP_release
```

where the folder names with `_FS` (e.g. `Sub0017_Ses1_FS`) meaning they are structual MRI data after `recon-all` processing, and the folder names without `_FS` (e.g. `Sub0017_Ses1`) meaning they are raw functional MRI data.

- Preprocessed data (groud truth)

4 GSP subjects were chosen for the Unit Tests and the data is stored here.

```
/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/FCmetrics
```

----

## Code

The users can use this command to call the preprocessing pipeline

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/FCmetrics/CBIG_preproc_unit_tests_FCmetrics_fMRI_preproc.csh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/FCmetrics/Preproc_Unit_Test_Pipeline.txt
```

The estimated walltime and mem usage are ~2h and ~2G.

----

## Results

The users can use the following MATLAB code to compare the FC metrices

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/FCmetrics/CBIG_preproc_FCmetrics_UnitTestComparison.m

```
The `CBIG_preproc_FCmetrics_UnitTestComparison.m` is a MATLAB function which takes in the <preproc_out_dir> as its input.

```

CBIG_preproc_FCmetrics_UnitTestComparison(<preproc_out_dir>)

```

If there are differences between any of the FC matrices with the ground truth, MATLAB will display a message "Difference noted in one or more correlation matrices. Check the inequal_corr_log.txt file in your testdir for more details". A text file will be created which shows the subject and which correlation matrix has an error. If there are no differences, MATLAB will display a message "All correlation matrices are the same."

```
