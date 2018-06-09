This readme includes the steps about how to run the single subject preprocessing pipeline to check the correctness of code (this check is especially useful for abmin when you are releasing some code). 

Notice that all the filenames and directories below work for CBIG lab only. It is assumed that the input surface data are in `fsaverage5` space.

----

## Data

- Data set

The `CBIG_fMRI_Preproc2016` preprocessing scripts use the Brain Genomics Superstruct Project (GSP) data as a baseline data set. All data mentioned here are based on GSP data set.

GSP data contain both structual MRI (T1) and functional MRI (T2*). All subjects (N = 1570) are healthy, young subjects (age: 18-35). The preprocessed (by `recon-all`) structual MRI data and raw functional MRI are in this folder:

```
/mnt/eql/yeo1/data/GSP_release
```

where the folder names with `_FS` (e.g. `Sub1116_Ses1_FS`) meaning they are structual MRI data after `recon-all` processing, and the folder names without `_FS` (e.g. `Sub1116_Ses1`) meaning they are raw functional MRI data.

- Preprocessed data (groud truth)

The preprocessed subject (Sub1116_Ses1) for comparison is stored here

```
/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject/data
```

----

## Code

The users can use this command to call the preprocessing pipeline

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_call_fMRI_preproc.csh <preproc_out_dir>
```

The pipeline will be run with this configuration file

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/prepro.config
```

The estimated walltime and mem usage are ~9h and ~30G.

----

## Results

The users can use 

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_cmp_in_fsaverage5.csh <preproc_out_dir> <file_stem_in_fsaverage5_space> <compare_out_dir>
```

to compare their results with the groud truth output mentioned in **Data** section in fsaverage space, and use

```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/unit_tests/single_subject/CBIG_preproc_unit_tests_cmp_in_mni2mm.csh <preproc_out_dir> <file_stem_in_fsaverage5_space> <compare_out_dir>
```

to compare their results with the groud truth output mentioned in **Data** section in MNI 2mm volumetric space.

----

## References

- Holmes, Avram J., et al. "Brain Genomics Superstruct Project initial data release with structural, functional, and behavioral measures." Scientific data 2 (2015).
