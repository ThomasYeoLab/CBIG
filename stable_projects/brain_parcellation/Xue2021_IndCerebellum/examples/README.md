# Examples of Individual Cerebellar Parcellation

----

References
==========

- Xue A, Kong R, Yang Q, et al. The Detailed Organization of the Human Cerebellum Estimated by Intrinsic Functional Connectivity Within the Individual. Journal of Neurophysiology, 2021

----

Data
====

In this example, we use the functional connectivity profiles of CoRR-HNU datasets on the cerebral cortical surface:
- **CoRR-HNU:**
  
  `$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0?/subj0?_sess?`

The cerebellar data for each subject is fake data in MNI 4mm space. 

----

Code
====

This example will generate 2 individual-specific cerebellar parcellation. The cerebral cortical parcellation are generated using [MS-HBM model](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Kong2019_MSHBM) based on a group-level cerebral cortical parcellation.

----

Run
====

- We provide a wrapper script which runs through the entire example:

`$CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/examples/CBIG_IndCBM_example_wrapper.m`

User can use this wrapper script as a reference. Please run the following command:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/examples
CBIG_IndCBM_example_wrapper(<output_dir>)
```
Results will be saved under `<output_dir>`.

The example includes a few steps:

### Create example data list

6 files will be genereated under `<output_dir>/list`:

```
lh_sub1_list.txt
lh_sub2_list.txt
rh_sub1_list.txt
rh_sub2_list.txt
vol_sub1_list.txt
vol_sub2_list.txt
```

Each subject should have 3 lists. 

`lh_sub?_list.txt` contains runs of left hemisphere time series (on surface).

`rh_sub?_list.txt` contains runs of right hemisphere time series (on surface).

`vol_sub?_list.txt` contains runs of cerebellum time series (in volume).

Each row of these files is the file path of a run. 

Files under `<output_dir>/MSHBM` follow the same structure with [MS-HBM model](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Kong2019_MSHBM). This is unnecessary if you are using your own individual-specific cereberal cortical parcellation.

### Cerebral cortical parcellation

This step first computes surface profile on fsaverage5 mesh and then average the profile between 2 subjects. Then average profiles are used to generate MS-HBM parameters for the given group-level parcellation `./group_surf/fsaverage5_10networks.mat`. Finally, use [MS-HBM model](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Kong2019_MSHBM) to generate individual-specific cereberal cortical parcellations. Results will be saved under `<output_dir>/MSHBM/ind_parcellation`

If you are using your own individual-specific cereberal cortical parcellation, make sure your parcellation of each subject includes 2 column vectors `lh_labels` and `rh_labels`. 

Please see `./examples/example_surf_labels/Sub1_example_surf_labels.mat` as a 10-network example on `fsaverage5`. 

### Create cifti template

Cifti template will be created based on fsaverge5 surface and volume masks `./input/mask/sub?/sub?_bin_mask_4mm.nii.gz`. 

The template file will be saved at `<output_dir>/template`. This template will be used in computing volume to surface functional connectivity and cerebellar parcellation.

### Cerebellar parcellation

Functional connectivity between cerebellum and cerebral cortex will be computed by `CBIG_IndCBM_compute_vol2surf_fc.m` and pass to the cerebellar parcellation function `CBIG_IndCBM_cerebellum_parcellation.m`

Final output of this step will be saved under `<output_dir>/parcellation`.

5 files will be genereated under `<output_dir>/parcellation/sub?`:

```
IndCBM_parcellation_top100.dlabel.nii
IndCBM_parcellation_top100_confidence.dscalar.nii
IndCBM_parcellation_top100.nii.gz
IndCBM_parcellation_top100_confidence.nii.gz
Freesurfer_cerebellum_LUT.txt
```

`IndCBM_parcellation_top100.dlabel.nii` is the final individual-specific parcellation including cerebral cortial surface and the cerebellum. 

`IndCBM_parcellation_top100_confidence.dscalar.nii` is the confidence of the cerebellar parcellation. Measured by `1 - N_second_frequent_network/N_most_frequent_network`

`IndCBM_parcellation_top100.nii.gz` is the cerebellar parcellation in nifti format. This file only covers the cerebellum..

`IndCBM_parcellation_top100_confidence.nii.gz` is the confidence of cerebellar parcellation in nifti format.

`Freesurfer_cerebellum_LUT.txt` is the freesurfer look up table of the cerebellar parcellation. Only used when visualizing `IndCBM_parcellation_top100.nii.gz` in freeview. 

- We also provide a script which checks user's example results with our reference results:

`$CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/examples/CBIG_IndCBM_check_example_results.m`

Assume your example results are saved in `<output_dir>`. Please run the following command to check your results:

```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/examples
CBIG_IndCBM_check_example_results(<output_dir>)
```

You will receive this message if your results are correct:
`Your example results are correct!`.

----

Bugs and Questions
====
Please contact XUE Aihuiping at xueaihuiping@gmail.com.

