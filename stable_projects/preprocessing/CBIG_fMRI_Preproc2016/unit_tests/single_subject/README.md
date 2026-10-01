This readme describes how to run the single subject preprocessing pipeline unit tests to check the correctness of code (this check is especially useful for admin when releasing code). Notice that all filenames and directories below work for CBIG lab only. It is assumed that the input surface data are in `fsaverage5` space.

In this unit test, there are 4 test cases:

1. Motion filtering (See [Respiratory pseudomotion filtering](https://github.com/YeoPrivateLab/CBIG_private/blob/develop/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/respiratory_pseudomotion_filtering.md) README)
2. aCompCor regression
3. Spatial distortion correction using opposite phase encoding direction fieldmaps
4. Motion correction (multi-echo, 2 runs, no STC/SDC)

The full preprocessing pipeline test (sub005, multi-echo, GSR regression) has been moved to the stepwise unit tests folder.

----

## Data
- Data set

For the **motion filtering**, **SDC opposite PED**, and **simple ME case** tests, the subject is from the Healthy Brain Network (**HBN**).

HBN data contain both structural MRI (T1) and functional MRI (T2*). All subjects (N = 2195) are young subjects (age: 6-22). The preprocessed (by `recon-all`) structural MRI data and raw functional MRI are in this folder:
`/mnt/nas/CSC7/Yeolab/Data/HBN`
where folder names with `_FS` (e.g. `sub-NDARBF851NH6_FS`) are structural MRI data after `recon-all` processing, and folder names without `_FS` (e.g. `sub-NDARBF851NH6`) are raw functional MRI data.

For the **aCompCor** test, the subject is from the Brain Genomics Superstruct Project (**GSP**).

GSP data contain both structural MRI (T1) and functional MRI (T2*). All subjects (N = 1570) are healthy, young subjects (age: 18-35). The preprocessed (by `recon-all`) structural MRI data and raw functional MRI are in this folder:
`/mnt/isilon/CSC2/Yeolab/Data/GSP`
where folder names with `_FS` (e.g. `Sub1116_Ses1_FS`) are structural MRI data after `recon-all` processing, and folder names without `_FS` (e.g. `Sub1116_Ses1`) are raw functional MRI data.


- Preprocessed data (ground truth)

  * For the **motion filtering** test, the preprocessed subject (`sub-NDARBF851NH6`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/test_motion_filter/sub-NDARBF851NH6`

  * For the **aCompCor** test, the preprocessed subject (`Sub1116_Ses1`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/aCompCor_deoblique/Sub1116_Ses1`

  * For the **SDC opposite PED** test, the preprocessed subject (`sub-NDARBF851NH6`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/sub-NDARBF851NH6`

  * For the **simple ME case** test, the preprocessed subject (`TMR00032_TMSD`) for comparison is stored here

  `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/TMR00032_TMSD`

----

## Code
The users can use the following commands to call each test case:

* Motion filtering:

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_test_motion_filtering.csh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro_motion_filtering.config
```

* aCompCor:

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc_aCompCor_deoblique.csh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro_aCompCor_deoblique.config
```

* SDC opposite PED:

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc_sdc_oppo_PED.csh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro_sdc_oppo_PED.config
```

* Motion correction:

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc_motion_correction.sh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro_motion_correction.config
```

**Alternatively**, users can run all 4 test cases in one wrapper function in MATLAB:

```
runtests('CBIG_preproc_single_subject_unit_test')
```

----

## Results

For the **motion filtering** test, only the motion parameters are compared with the ground truth:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/test_motion_filter/sub-NDARBF851NH6/bold/<run_number>/sub-NDARBF851NH6_bld<run_number>_rest_skip8_stc_mc.par`

For the **aCompCor** test, the residualised volumes are compared with the ground truth:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/aCompCor_deoblique/Sub1116_Ses1/bold/<run_number>/Sub1116_Ses1_bld<run_number>_rest_skip4_stc_mc_resid.nii.gz`

For the **SDC opposite PED** test, the distortion-corrected volumes are compared with the ground truth:

`$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/sub-NDARBF851NH6/bold/<run_number>/sub-NDARBF851NH6_bld<run_number>_rest_skip8_stc_mc_sdc.nii.gz`

For the **simple ME case** test, the MC volumes, final bandpass volumes, and surface time series are compared with the ground truth in `$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data/TMR00032_TMSD`. All comparisons use a minimum voxel/vertex correlation threshold of 0.99.

----

## References
- Holmes, Avram J., et al. "Brain Genomics Superstruct Project initial data release with structural, functional, and behavioral measures." Scientific data 2 (2015).
- Alexander, L. et al. An open resource for transdiagnostic research in pediatric mental health and learning disorders. Scientific Data 4, 170181 (2017).
