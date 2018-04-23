This README includes the instruction on how to run CBIG_RF unit test. Notice that all filenames and directories in this unit test **work for CBIG lab only**

----

## Reference

Jianxiao Wu, Gia H. Ngo, Alexander Schaefer, Douglas Greve, Bruce Fischl, Simon B. Eickhoff, B.T. Thomas Yeo. **Accurate Nonlinear Mapping between MNI152/Colin27 Volumetric and FreeSurfer Surface Coordinate Systems**. Under review.

----

## Data

The unit test script uses the Brain Genomics Superstruct Project (GSP) data. The preprocessed (by `recon-all`) structural MRI data are in the folder:
```
/mnt/eql/yeo1/data/GSP_release
```
where each subject folder is named by subject id + `_FS` (e.g. `Sub0001_Ses1_FS`)

In total, 5 subjects were used in the unit test, corresponding to the first 5 subjects using the default subject id list in `CBIG_RF` scripts (`$CBIG_CODE_DIR/stable_projects/registration/Wu2017_RegistrationFusion/bin/GSP_subjectid.csv`).

In addition, the `warp/` folder contains the ANTs registration warps between those 5 subjects and MNI152 space.

----

## Code

The codes for unit test are in the folder:

`/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/registration/Wu2017_RegistrationFusion`

- MNI152-to-fsaverage test

To run the unit test, call `CBIG_RF_unit_test_vol2surf.sh`, using the `-o` option to specify where the results should be put.

The final projection results (in left hemisphere) will be compared to the default results in `data/` folder. The unit test is successful if the screen prints `The two volumes are identical`.

You can also check the right hemisphere results by comparing your projected result with `data/rh.projected_central_sulc.nii.gz`

- fsaverage-to-MNI152 test

To run the unit test, call `CBIG_RF_unit_test_surf2vol.sh`, using the `-o` option to specify where the results should be put.

The final projection results will be compared to the default results in `data/` folder. The unit test is successful if the screen prints `The two volumes are identical`.



