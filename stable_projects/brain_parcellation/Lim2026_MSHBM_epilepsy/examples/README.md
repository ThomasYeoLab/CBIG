# Examples: Individual-Specific Epilepsy Network Parcellation

This folder demonstrates how to generate an individual-specific 15-network cortical
parcellation for a drug-resistant epilepsy patient from resting-state fMRI data, using
the NIH 34-subject MS-HBM epilepsy model. The script will also compute the **language laterality
index (LI)** for the individual subject.

---

## Quick Start

In MATLAB, run:

```matlab
CBIG_MSHBM_Epilepsy_wrapper();
```

This will process the 5 example fMRI runs in `input/` and save the parcellation, a
brain surface visualisation, and the language laterality index to `out_dir/`.

To check that the outputs match the reference results:

```matlab
CBIG_MSHBM_Epilepsy_check_example_results('out_dir');
```

The reference outputs are stored in `results/` and should not be overwritten. Always
pass a separate output directory to the wrapper and to the check function.

---

## Folder Contents

```
examples/
├── CBIG_MSHBM_Epilepsy_wrapper.m                    ← Example script using provided input data
├── CBIG_MSHBM_Epilepsy_check_example_results.m      ← Checks output against reference
├── input/
│   ├── lh.example_bld001_fs6_sm6.nii.gz             ← Left hemisphere fMRI, run 1
│   ├── lh.example_bld002_fs6_sm6.nii.gz             ← Left hemisphere fMRI, run 2
│   ├── lh.example_bld003_fs6_sm6.nii.gz             ← (runs 3–5 similarly named)
│   ├── rh.example_bld001_fs6_sm6.nii.gz             ← Right hemisphere fMRI, run 1
│   ├── ...
│   ├── example_bld001_motion_outliers.txt            ← Motion censor file, run 1
│   └── ...                                           ← (runs 2–5 similarly named)
├── out_dir/                                          ← Created on first run (default output)
│   ├── data_list/
│   │   ├── censor_list/sub1_sess{1-5}.txt            ← Censor file lists per session
│   │   └── fMRI_list/{lh,rh}_sub1_sess{1-5}.txt      ← fMRI file lists per session
│   ├── profile_list/test_set/{lh,rh}_sess{1-5}.txt   ← FC profile path lists
│   ├── profiles/sub1/sess{1-5}/
│   │   ├── lh.sub1_sess{N}_fsaverage6_roifsaverage3.surf2surf_profile.nii.gz
│   │   └── rh.sub1_sess{N}_fsaverage6_roifsaverage3.surf2surf_profile.nii.gz
│   ├── priors/                                       ← Symlink to group prior directory
│   ├── ind_parcellation/test_set/
│   │   └── Ind_parcellation_MSHBM_sub1_w80_MRF10.mat ← Individual parcellation labels
│   ├── epilepsy_w80c10.png                           ← Brain surface visualisation
│   └── LI.mat                                        ← Language laterality index (LI_lang)
└── results/                                          ← Reference outputs (do not overwrite)
    ├── ind_parcellation/test_set/
    │   └── Ind_parcellation_MSHBM_sub1_w80_MRF10.mat
    ├── epilepsy_w80c10.png
    └── LI.mat
```

The core parcellation and LI computation are in `CBIG_MSHBM_Epilepsy_LI.m` at the
top level of this release. The wrapper script is a thin example that sets up the
input paths and calls it.

---

## Key Parameters

**`target_mesh`** — surface mesh resolution of your fMRI data

| Option | Description |
|---|---|
| `'fsaverage6'` | FreeSurfer fsaverage6 |

The NIH 34-subject MS-HBM epilepsy model was trained on fsaverage6 data. Using any
other mesh resolution requires retraining the group prior and is not supported by
this release.

**`censor_list`** — motion censoring specification

| Option | Description |
|---|---|
| Cell array of file paths | One `.txt` file per run; each file has one row per timepoint (1 = keep, 0 = censor) |
| `'NONE'` | No censoring — use when fMRI files are already frame-filtered |

**`w`** — group spatial prior weight

Higher values pull the individual parcellation closer to the group-average network
topology. The NIH 34-subject MS-HBM epilepsy model uses `w = 80`.

**`c`** — MRF smoothness constraint weight

Higher values produce spatially smoother parcellations. The NIH 34-subject MS-HBM epilepsy model uses `c = 10`.

---

Please contact Mervyn Lim Jun Rui at [mervynlim@u.nus.edu](mailto:mervynlim@u.nus.edu) for any questions.
