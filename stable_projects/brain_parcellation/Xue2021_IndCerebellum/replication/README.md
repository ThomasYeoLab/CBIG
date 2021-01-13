# Replication of Individual Cerebellar Parcellation

----

Notice that all filenames and directories in this replication **work for CBIG lab only**.

References
==========
- Xue A, Kong R, Yang Q, et al. The Detailed Organization of the Human Cerebellum Estimated by Intrinsic Functional Connectivity Within the Individual. Journal of Neurophysiology, 2021

----

Data
====

2 individual subjects scanned at the Harvard Center for Brain Science. 

31 sessions collected for each subject. Data is divided into discovery set and replication set by odd and even sessions. 

- Discovery set: 16 sessions
- Replication set: 15 sessions
- All set: 31 sessions

----

Code
====

The replication wrapper `CBIG_IndCBM_replication_wrapper.sh` will generate individual-specific cerebellar parcellation for each subject each set. 

----

Run
====

First, please add `CBIG_IndCBM_REP_DIR` to your config file, refering to `./config/CBIG_IndCBM_tested_config.sh`.

Use the following wrapper script to runs through the entire replication:

`$CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/replication/CBIG_IndCBM_replication_wrapper.sh`

Note that the whole wrapper take ~20 hours to run and the output folder will take ~290GB space. 

Please run the following command:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/replication
sh CBIG_IndCBM_replication_wrapper.sh <output_dir>
```
Results will be saved under `<output_dir>`.

The replication includes a few steps:

### Create replication data list

5 files will be genereated under `<output_dir>/<subject>/<set>`:

```
sess_<set>_list.txt
runs_<set>_list.txt
lh_<set>_list.txt
rh_<set>_list.txt
MNI_<set>_list.txt
```

`sess_<set>_list.txt` contains session name for the corresponding set.

`run_<set>_list.txt` contains run name for the corresponding set.

Each row of these files is only the session/run name.

`lh_<set>_list.txt` contains runs of left hemisphere time series (on surface).

`rh_<set>_list.txt` contains runs of right hemisphere time series (on surface).

`MNI_<set>_list.txt` contains runs of cerebellum time series (in volume).

Each row of these files is the file path of a run. 

Files under `<output_dir>/<set>/MSHBM` follow the same structure with [MS-HBM model](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Kong2019_MSHBM). 

### Cerebral cortical parcellation

This step first computes surface profile on `fsaverage6` mesh and then average the profile between 2 subjects. Then average profiles are used to generate MS-HBM parameters for the given group-level parcellation `./input/HCP_10networks.mat`. Finally, use [MS-HBM model](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Kong2019_MSHBM) to generate individual-specific cereberal cortical parcellations. Results will be saved under `<output_dir>/<set>/MSHBM/ind_parcellation`

### Create cifti template

Cifti template will be created based on fsaverge6 surface and volume masks `$CBIG_IndCBM_REP_DIR/Derivatives/<sub>/mask/<sub>_bin_mask.nii.gz`. 

The template file will be saved at `<output_dir>/template`. This template will be used in computing volume to surface functional connectivity and cerebellar parcellation.

### Compute functional connectivity

Functional connectivity between cerebellum and cerebral cortex will be computed by `CBIG_IndCBM_compute_vol2surf_fc.m` and saved under `<output_dir>/FC`.

Note that each `.mat` file takes 45~50GB. So please pay attention to the remaining disk space.

### Cerebellar parcellation

The cerebellar parcellation function `CBIG_IndCBM_cerebellum_parcellation.m` will generate cerebellar parcellation based on the functional connectivity matrix in the above step.

Final output will be saved under `<output_dir>/parcellation`.

For each subject, each set, there will be 2 files:

```
<subject>_<set>_IndCBM_parcellation.dlabel.nii
<subject>_<set>_IndCBM_parcellation_confidence.dscalar.nii
```

`<subject>_<set>_IndCBM_parcellation.dlabel.nii` is the final individual-specific parcellation including cerebral cortial surface and the cerebellum. 

`<subject>_<set>_IndCBM_parcellation_confidence.dscalar.nii` is the confidence of the cerebellar parcellation. Measured by `1 - N_second_frequent_network/N_most_frequent_network`

Assume your replication results are saved in `<output_dir>`. Please run the following command to check your results:

```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/replication
CBIG_IndCBM_check_replication_results(<output_dir>)
```

You will receive this message if your results are correct:
`Your replication results are correct!`.

Reference files are saved at `$CBIG_IndCBM_REP_DIR/ref_replication_results`.

Job log files are saved under `<output_dir>/logs`.

----

Bugs and Questions
====
Please contact XUE Aihuiping at xueaihuiping@gmail.com.

