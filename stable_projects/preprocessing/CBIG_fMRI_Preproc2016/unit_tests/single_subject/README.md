This readme includes the steps about how to run the single subject preprocessing pipeline to check the correctness of code (this check is especially useful for admin when you are releasing some code). Notice that all the filenames and directories below work for CBIG lab only. It is assumed that the input surface data are in `fsaverage5` space.

In this unit tests, there are in total 3 different test cases:

1. Full preprocessing pipeline test (regression: GSR; spatial distortion correction: opposite phase-encoding direction)
2. aCompCor test (regression: aCompCor)
3. spatial distortion correction using magnitude and phase difference test (spatial distortion correction: magnitude and phase difference)

----

## Data

- Data set

In this unit test, we are using 3 subjects from 2 datasets. 

For test case **1** and **3**, these 2 subjects are from the Healthy Brain Network (**HBN**).

HBN data contain both structual MRI (T1) and functional MRI (T2*). All subjects (N = 2195) are young subjects (age: 6-22). The preprocessed (by `recon-all`) structual MRI data and raw functional MRI are in this folder:

```
/mnt/eql/yeo11/data/HBN
```

where the folder names with `_FS` (e.g. `sub-NDARBF851NH6_FS`) meaning they are structual MRI data after `recon-all` processing, and the folder names without `_FS` (e.g. `sub-NDARBF851NH6`) meaning they are raw functional MRI data.



For test case **2**, the subject is from the Brain Genomics Superstruct Project (**GSP**).

GSP data contain both structual MRI (T1) and functional MRI (T2*). All subjects (N = 1570) are healthy, young subjects (age: 18-35). The preprocessed (by `recon-all`) structual MRI data and raw functional MRI are in this folder:

```
/mnt/eql/yeo1/data/GSP_release
```

where the folder names with `_FS` (e.g. `Sub1116_Ses1_FS`) meaning they are structual MRI data after `recon-all` processing, and the folder names without `_FS` (e.g. `Sub1116_Ses1`) meaning they are raw functional MRI data.



- Preprocessed data (groud truth)

  * For test case 1 (full preprocessing pipeline), the preprocessed subject (`sub-NDARBF851NH6`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/sub-NDARBF851NH6`

  * For test case 2 (aCompCor), the preprocessed subject (`Sub1116_Ses1`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/aCompCor/Sub1116_Ses1`

  - For test case 3 (spatial distortion correction magnitude and phase difference method), the preprocessed subject (`sub-NDARAA075AMK`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/sub-NDARAA075AMK`

----

## Code

The users can use the following commands to call each of the 3 test cases:

* Test case 1 (full preprocessing pipeline):

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc.csh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro.config
```

* Test case 2 (aCompCor):

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc_aCompCor.csh
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro_aCompCor.config
```

* Test case 3 (spatial distortion corection magnitude and phase difference method):

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc_sdc_mag_phasediff.csh
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro_sdc_mag_phasediff.config
```



**Alternatively**, users can run these 3 test cases in one wrapper function as following in MATLAB.

```
runtests('CBIG_preproc_single_subject_unit_test')
```

The estimated walltime and mem usage are ~1.5h and ~6G.

----

## Results

For test case 1 (full preprocessing pipeline), the users should compare **three** aspects of their results from the ground truth: (1) the timeseries in fsaverage surface space; (2) the timeseries in MNI 2mm volumetric space; (3) the functional connectivity matrices

- Compare timeseries in fsaverage space

In a bash shell, use the following command

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_cmp_in_fsaverage5.csh <preproc_out_dir> <file_stem_in_fsaverage5_space> <compare_out_dir>
```

This command will generate the correlation between the timeseries in your preprocessed NIFTI file and the ground truth timeseries for each vertex in fsaverage5 space. You can check the histogram of the correlations (named as `<compare_out_dir>/sub-NDARBF851NH6/<run_number>/gt_user-test_corr_surf_hist.png`), as well as the plot of correlations projected on fsaverage5 space (named as `<compare_out_dir>/sub-NDARBF851NH6/<run_number>/gt_user-test_corr_surf.png`).

- Compare timeseries in MNI 2mm space

In a bash shell, use the following command:

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_cmp_in_mni2mm.csh <preproc_out_dir> <final_file_stem_in_MNI_space> <compare_out_dir>
```

This command will generate the correlation between the timeseries in your preprocessed NIFTI file and the ground truth timeseries for each grey matter voxel in MNI 2mm space. You can check the histogram of the correlations (named as `<compare_out_dir>/sub-NDARBF851NH6/<run_number>/gt_user-test_corr_vol_gm_hist.png`), as well as visualize the spatial distribution of the correlations in the MNI space (named as `<compare_out_dir>/sub-NDARBF851NH6/<run_number>/gt_user-test_corr_vol_gm.nii.gz`) using freeview.

- Compare functional connectivity matrices

In a Matlab terminal, use the following command:

```
cd(fullfile(getenv('$CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'));
CBIG_preproc_FCmatrices_UnitTestCmp('<preproc_out_dir>', '<test_file_stem>');
```

`preproc_out_dir` should store the directory path where the output data (generated by `CBIG_preproc_unit_tests_call_fMRI_preproc.csh`) is stored. `test_file_stem` is an optional argument. It is a string that is appended to the generated FC matrices filename. For example, a typical FC matrix file name could be `sub-NDARBF851NH6_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_all2all.mat`. In this case, the `test_file_stem` would be `rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6`. The default `test_file_stem` is `rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6`.  If you successfully replicated the ground truth functional connectivity matrices, this command will output `"All correlation matrices are the same"`. A text file will named `inequal_corr_log.txt` will also be generated and stored in `preproc_out_dir`. It will display the differences between the user-generated results and the ground truth.



For test case 2 and 3 (aCompCor & spatial distortion correction magnitude and phase difference method), only the final volume is compared with the ground truth. 

Specifically, for test case 2, the volumes that are checked include: `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/aCompCor/Sub1116_Ses1/bold/<run_number>/Sub1116_Ses1_bld<run_number>_rest_skip4_stc_mc_resid.nii.gz` 

For test case 3, the volume that is checked include:`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/sub-NDARAA075AMK/bold/001/sub-NDARAA075AMK_bld001_rest_skip8_stc_mc_sdc.nii.gz`

----

## References

- Holmes, Avram J., et al. "Brain Genomics Superstruct Project initial data release with structural, functional, and behavioral measures." Scientific data 2 (2015).
- Alexander, L. et al. An open resource for transdiagnostic research in pediatric mental health and learning disorders. Scientific Data 4, 170181 (2017).
