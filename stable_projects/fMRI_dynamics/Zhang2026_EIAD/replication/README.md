## Replication

This folder contains the scripts and participant tables needed to reproduce the
main results of the paper in the ADNI and MACC cohorts.

### Scripts

| Script | Purpose |
| --- | --- |
| `CBIG_EIAD_ADNI.m` | All main ADNI analyses. Writes `../output/ADNI.mat`. |
| `CBIG_EIAD_MACC.m` | All main MACC analyses. Writes `../output/MACC.mat`. |
| `CBIG_EIAD_run_harmonization.m` | ComBat harmonization of the regional E/I ratios across scanners. |
| `CBIG_EIAD_draw_surface_map.m` | Renders a per-region value map on the fsaverage5 surface. |

`CBIG_EIAD_ADNI.m` and `CBIG_EIAD_MACC.m` each run three analyses:

1. a GLM between the mean cortical E/I ratio and the amyloid marker (CSF Abeta42
   in ADNI, plasma p-tau217 in MACC), adjusting for age, sex, head motion and
   years of education;
2. the same GLM region by region, and the Spearman correlation of the resulting
   68 slopes with the sensorimotor-association (SA) axis;
3. a mediation analysis of the amyloid marker -> cognition relationship with the
   E/I ratio as the mediator, at the mean-cortex and regional levels.

Cognition is the ADSP-PHC memory composite in ADNI and a single-factor score over
the verbal-memory and MMSE items in MACC.

### Usage

> **The released tables are not directly runnable.** To comply with the ADNI and
> MACC data-use agreements, `ADNI_df.csv` and `MACC_df.csv` ship with only the
> subject IDs and the model-derived neural columns (the regional E/I ratios, plus
> the underlying S_E and S_I) populated. Every demographic, biomarker, and
> cognitive column is present as a header but left empty. Before the wrappers will
> run, fill these columns in from your own ADNI / MACC download using the released
> subject IDs — see [Repopulating the tables](#repopulating-the-tables) below.

From the `script` folder, run

```matlab
CBIG_EIAD_ADNI
CBIG_EIAD_MACC
```

in MATLAB. Results are written to `../output/` and can be compared against
`../reference_output/`. The regional mediation dominates the run time (68 regions
x 10,000 bootstrap draws per cohort); each wrapper takes roughly 2 minutes on
a desktop CPU. To check your setup quickly, run the example in `../../examples/`
instead, which takes under 10 seconds.

These scripts start from the participant tables released here, which already
contain the harmonized regional E/I ratios. Fitting the pFIC model and running
the harmonization are separate, earlier steps and are documented below.

Harmonization is run separately, since it operates on the participant-level
simulation output rather than the tables released here:

```matlab
EI_harmonized = CBIG_EIAD_run_harmonization(scanner_list_dir, cov_file, ...
    S_E_path, S_I_path, r_E_path, combat_dir);
```

See the function header for the meaning of each argument. `S_E_path`, `S_I_path`
and `r_E_path` are path templates containing a `%s` placeholder for the subject
ID.

### From pFIC output to regional E/I ratio

`CBIG_EIAD_run_harmonization.m` implements the full path from the pFIC
simulation output to the `HarmonizedEI*` columns of `ADNI_df.csv`, in two
steps.

**1. E/I ratio.** For each participant, the pFIC fit provides the regional
excitatory (`S_E`) and inhibitory (`S_I`) synaptic gating variables as
`68 x #simulation` matrices — one value per Desikan region per simulation.
Simulations whose excitatory firing rate `r_E` falls outside the physiological
window (2.7–3.3 Hz) are discarded, and the surviving simulations are averaged to
give one `S_E` and one `S_I` value per region.

The regional E/I ratio is then the elementwise ratio of these two vectors,
`S_E ./ S_I`.

**2. ComBat.** The resulting 68 x N matrix is passed to `combat.m` with scanner
identity as the batch variable. The biological covariates retained in the model
matrix are age, sex, diagnosis (MCI and AD indicators, CN as reference), the two
age-by-diagnosis interactions, and mean framewise displacement, so that ComBat
removes scanner effects while preserving variance attributable to these
variables. The non-parametric empirical Bayes adjustment is used (the trailing
`0` argument to `combat`). MACC is single-site and is therefore not harmonized;
`MACC_df.csv` holds the unharmonized ratios in `EI1` ... `EI68`.

Scanner assignment comes from the `.txt` files in `scanner_list_dir`, one file
per scanner. Files are read in the order returned by `dir()`, and the k-th file
defines batch k; the batch labels are arbitrary, so this ordering does not
affect the result.

### Mediation analysis

The mediation is a nonparametric bootstrap, implemented in the local function
`bootstrap_mediation` at the bottom of `CBIG_EIAD_ADNI.m`; the same procedure is
used in `CBIG_EIAD_MACC.m` and in the example, so that each script runs standalone.

The independent variable X is the amyloid marker (log10 CSF Abeta42 in ADNI,
log10 plasma p-tau217 in MACC), the mediator M is the E/I ratio (mean cortical, or a
single region), and the dependent variable Y is the memory composite. X, M and Y
are z-scored, so all paths are standardized coefficients. Age, sex, mean
framewise displacement and years of education enter every model as covariates,
with the continuous ones z-scored.

Three ordinary least-squares models are fitted on the original sample:

| Model | Coefficient of interest |
| --- | --- |
| `M ~ X + Z` | path a |
| `Y ~ M + X + Z` | path b, and the direct effect c' (the X coefficient) |
| `Y ~ X + Z` | total effect c |

The indirect effect is the product `a * b`.

All three models are then refitted on 10,000 bootstrap samples drawn with
replacement at the participant level (`randsample(n, n, true)`), giving a
bootstrap distribution for each path. The p-value for a path is two-sided and
read off that distribution as `2 * min(P(boot >= 0), P(boot <= 0))`; the
confidence interval for the indirect effect is the 2.5th-97.5th percentile
range. The seed is fixed at the top of the function, so the reported
values reproduce exactly.

The regional analysis repeats this for each of the 68 regions in turn, and the
resulting p-values are corrected across regions with Benjamini-Hochberg FDR
(`mafdr(..., 'BHFDR', 1)`); the calls are in the commented-out visualization
block. This regional loop (68 regions x 10,000 draws) dominates the run time of
both wrappers.

### Data files

`ADNI_df.csv` (302 participants) and `MACC_df.csv` (240 participants) hold one
row per participant.

To comply with the ADNI and MACC data-use agreements, these tables release real
values for only two things: the **subject ID** and the **model-derived neural
columns** (the regional E/I ratios, plus the underlying regional S_E and S_I).
The `ID` column holds the original ADNI / MACC subject identifiers, so each row
can be linked back to the source cohort by anyone with approved access. Every
other column — all demographics, biomarkers, and cognitive scores — is kept as a
header but left **empty**; those values are not redistributed here.

The column tables below document what each field means and mark whether it is
**released** or shipped **empty (fill from ADNI / MACC)**. See
[Repopulating the tables](#repopulating-the-tables) for how to fill the empty
columns so the wrappers run.

Columns shared by both tables:

| Column | Description | Released? |
| --- | --- | --- |
| `ID` | original subject ID (links to the source cohort) | **released** |
| `age` | age at the fMRI session (years) | empty |
| `sex_F0M1` | 0 = female, 1 = male | empty |
| `FD` | mean framewise displacement (mm) | empty |
| `yrsOfEducation` | years of education | empty |
| `BMI`, `Hypertension`, `Diabetes`, `Smoker` | vascular risk factors | empty |

ADNI-specific columns:

| Column | Description | Released? |
| --- | --- | --- |
| `AB42_CSF`, `AB40_CSF` | CSF Abeta42 and Abeta40 (pg/mL) | empty |
| `amyloid_status` | amyloid positivity | empty |
| `PHCMemory` | ADSP-PHC memory composite | empty |
| `composite_memory_score` | memory composite matched to the MACC items | empty |
| `HarmonizedEI1` ... `HarmonizedEI68` | ComBat-harmonized regional E/I ratios | **released** |
| `UnharmonizedEI1` ... `UnharmonizedEI68` | the same ratios before harmonization | **released** |
| `HarmonizedSE1` ... `HarmonizedSE68` | harmonized regional S_E | **released** |
| `HarmonizedSI1` ... `HarmonizedSI68` | harmonized regional S_I | **released** |

MACC-specific columns:

| Column | Description | Released? |
| --- | --- | --- |
| `ptau217`, `ptau181` | plasma p-tau217 and p-tau181 (pg/mL) | empty |
| `vbmwlrimdtrecal`, `vbmwlrdelyrecal` | word-list immediate and delayed recall | empty |
| `vbmsraimdt`, `vbmsradelyrecal` | story immediate and delayed recall | empty |
| `orientday_BL` ... `orientcountry_BL` | the 10 MMSE orientation items | empty |
| `registrationscore_BL`, `recall_BL` | MMSE registration and recall | empty |
| `EI1` ... `EI68` | regional E/I ratios | **released** |
| `SE1` ... `SE68`, `SI1` ... `SI68` | the underlying regional S_E and S_I | **released** |

MACC is single-site, so the MACC E/I ratios are not harmonized. Regions follow
the 68-region Desikan parcellation; `SA_rank_desikan.txt` gives the SA-axis rank
of each region in the same order.

### Repopulating the tables

The wrappers read every column by name, so to reproduce the published results you
only need to fill the empty columns back in — the row order, the `ID` column, and
all of the E/I / S_E / S_I columns are already in place and should not be changed.

With approved access to the source data:

1. Obtain the raw variables from the source cohort. ADNI data are available from
   https://adni.loni.usc.edu under the ADNI Data Use Agreement; MACC data are
   available from the Memory Ageing and Cognition Centre on reasonable request.
   The Released? tables above name every field you need, together with what it
   is (e.g. `AB42_CSF` = CSF Abeta42 in pg/mL, `PHCMemory` = ADSP-PHC memory
   composite, `ptau217` = plasma p-tau217, and so on).
2. Match each released `ID` to the corresponding participant in your download and
   copy the values into the matching empty columns, leaving the header names and
   the row order unchanged. Do not reorder or drop rows: the E/I / S_E / S_I
   columns are aligned to the existing row order.
3. Keep the units and coding shown in the column tables (for example
   `sex_F0M1` = 0 for female and 1 for male, and CSF / plasma markers in pg/mL),
   since the scripts z-score and `log10`-transform specific columns and assume
   these conventions.

Some ADNI / MACC variables have more than one record per participant (for
example, repeated CSF draws or cognitive assessments across visits). In that
case the value used in the paper is the one closest in time to the participant's
fMRI session, and picking a different visit will change the results. If you are
unsure which record to use, please contact the corresponding author
(Shaoshi Zhang, 0zhangshaoshi0@gmail.com; Thomas Yeo yeoyeo02@gmail.com) for 
the exact CSF, fMRI, and assessment dates used for each participant, so you can 
select the matching entry.

Once the columns are filled, run `CBIG_EIAD_ADNI` and `CBIG_EIAD_MACC` as
described under [Usage](#usage). The scripts themselves are unchanged from the
version used for the paper, so a correctly repopulated table reproduces the
values in `../reference_output/`.

### pFIC configuration files

`pFIC_ADNI_example_config.ini` and `pFIC_MACC_example_config.ini` are the configuration
files used to estimate the regional model parameters with the pFIC model
(https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/fMRI_dynamics/Zhang2024_pFIC/model). 
They record the cohort-specific acquisition settings (TR, number of frames, simulation period, 
sliding-window length) alongside the neural-mass, hemodynamic and training parameters.

The model was fit separately for each participant. Because the fit is
computationally intensive - approximately 5 hours per participant on a single
NVIDIA RTX 3090 GPU - it was run with a single noise instantiation per
participant. The random seed is recorded in the configuration files, so the fits can be
reproduced exactly.

The estimation procedure itself is not reimplemented here. The model code, the
training/validation/test workflow and the meaning of every field in these
configuration files are documented in the pFIC release
(https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/fMRI_dynamics/Zhang2024_pFIC),
which these files are intended to be used with.

### Requirements

MATLAB with the Statistics and Machine Learning Toolbox (`fitglm`, `factoran`,
`zscore`, `randsample`, `dummyvar`).

`CBIG_EIAD_run_harmonization.m` additionally needs `combat.m` from
https://github.com/Jfortin1/ComBatHarmonization (`Matlab/scripts`).

`CBIG_EIAD_draw_surface_map.m` additionally needs the CBIG repository
(https://github.com/ThomasYeoLab/CBIG) on the MATLAB path, together with
`1000subjects_clusters007_ref.mat` in the `script` folder. The calls to it are commented
out in both wrappers, so the analyses run without it.

The FDR-thresholded surface figures in both wrappers use `mafdr` (`'BHFDR', 1`)
from the Bioinformatics Toolbox. These calls sit in the commented-out
visualization blocks, so the Bioinformatics Toolbox is only required if you
uncomment them.

### Data availability

The demographic, biomarker, and cognitive values in `ADNI_df.csv` and
`MACC_df.csv` are not redistributed here; those columns ship empty and must be
repopulated from the source cohorts as described in
[Repopulating the tables](#repopulating-the-tables). Likewise, the
participant-level pFIC simulation output (`S_E`, `S_I`, `r_E`) and the covariate
and scanner lists used for harmonization are not redistributed. ADNI data are
available from https://adni.loni.usc.edu subject to the ADNI Data Use Agreement;
MACC data are available from the Memory Ageing and Cognition Centre on reasonable
request.
