# Replication

This folder contains scripts to replicate all figures and statistics from Lim 2026.

> **Note:** Replication scripts require access to the patient fMRI data, which is not
> included in this release. Input data may be obtained upon reasonable request to the authors.

---

## Environment Setup

Before running any replication scripts, configure your shell environment:

```bash
source replication/config/CBIG_epilepsy_config.sh
```

Edit `CBIG_epilepsy_config.sh` to set the paths to your CBIG installation, FreeSurfer,
MATLAB, SPM, AFNI, ANTs, Workbench, and FSL. The `CBIG_epilepsy_startup.m` file is loaded
automatically at MATLAB startup via `$MATLABPATH`.

---

## Folder Structure

```
replication/
├── config/
│   ├── CBIG_epilepsy_config.sh          Shell environment variables (edit before use)
│   └── CBIG_epilepsy_startup.m          MATLAB startup configuration
├── NIH/
│   ├── CBIG_train_MSHBM_Epilepsy.m        Step 1: Train MS-HBM (sequential; run directly or via submit)
│   ├── CBIG_train_MSHBM_Epilepsy_submit.m Step 1 (PBS): Submits the train script as a single PBS job
│   ├── CBIG_test_MSHBM_Epilepsy.m         Step 2: Test on 14 NIH subjects (Fig. 3A)
│   ├── CBIG_language_prediction.m        Step 3: Language LI prediction (Fig. 5)
│   ├── CBIG_ParcellationHomogeneity_FS_meantimecourse.m  Helper: FC homogeneity
│   ├── CBIG_supp_fig1_fc_similarity.m   Supp. Fig. 1: NIH–HCP FC similarity map
│   ├── CBIG_supp_fig3_dice.m            Supp. Fig. 3: Inter-subject DICE
│   ├── CBIG_supp_fig4_lang_homo.m       Supp. Fig. 4: Language network homogeneity
│   └── CBIG_supp_fig5_model_compare.m   Supp. Fig. 5: Homogeneity vs training set size
├── esfmri/
│   ├── CBIG_test_MSHBM_esfmri_homo.m   esfMRI resting-state homogeneity
│   ├── CBIG_test_MSHBM_esfmri_glm.m    esfMRI stimulation GLM gamma maps
│   └── CBIG_test_MSHBM_esfmri_inhomo.m esfMRI stimulation inhomogeneity
└── mshbm/
    ├── HCP/                              HCP prior and group parcellation
    ├── Du/                               Du prior and group parcellation
    ├── NIH/                              NIH prior and group parcellation
    └── union_medial_wall.mat             Medial wall mask for fsaverage6
```

Input fMRI data (not included in this release) is read from the path set by
`$CBIG_EPILEPSY_DATA_DIR` (defined in `config/CBIG_epilepsy_config.sh`):

```
$CBIG_EPILEPSY_DATA_DIR/
├── NIH/input/{subject}/surf/     — NIH preprocessed surface files
├── NIH/input/{subject}/qc/       — motion outlier files
├── NIH/input/FC_compare/         — pre-computed FC blocks (Supp. Fig. 1)
├── esfmri/input/rs_pp/{subject}/ — esfMRI resting-state surface files
└── esfmri/input/es_pp/{subject}/ — esfMRI stimulation surface files
```

Results are written to `replication/NIH/results/` and `replication/esfmri/results/`
within the CBIG codebase.

---

## Run Order

Each script is self-contained (sets up its own MATLAB path) and can be run from any
working directory in MATLAB.

### NIH cohort

| Order | Script | Output |
|-------|--------|--------|
| 1 | `NIH/CBIG_train_MSHBM_Epilepsy_submit.m` | Submits single PBS job; output in `NIH/results/mshbm_epilepsy/` |
| 2 | `NIH/CBIG_test_MSHBM_Epilepsy.m` | Fig. 3A: leave-one-run-out homogeneity |
| 3 | `NIH/CBIG_language_prediction.m` | Fig. 5: language LI boxplots and ROC curves |
| 4 | `NIH/CBIG_supp_fig1_fc_similarity.m` | Supp. Fig. 1: NIH–HCP FC similarity |
| 5 | `NIH/CBIG_supp_fig3_dice.m` | Supp. Fig. 3: inter-subject DICE |
| 6* | `NIH/CBIG_supp_fig4_lang_homo.m` | Supp. Fig. 4: language network homogeneity |
| 7* | `NIH/CBIG_supp_fig5_model_compare.m` | Supp. Fig. 5: model comparison |

\* Steps 4–7 require results from step 2 and the esfMRI scripts below.

### esfMRI cohort

| Order | Script | Output |
|-------|--------|--------|
| 1 | `esfmri/CBIG_test_MSHBM_esfmri_homo.m` | esfMRI per-subject homogeneity |
| 2 | `esfmri/CBIG_test_MSHBM_esfmri_glm.m` | Stimulation GLM gamma maps |
| 3 | `esfmri/CBIG_test_MSHBM_esfmri_inhomo.m` | Stimulation network inhomogeneity |

The esfMRI scripts can be run independently of the NIH scripts (they do not depend on
step 1 of the NIH cohort). Run esfMRI steps 1–3 before NIH step 6.
