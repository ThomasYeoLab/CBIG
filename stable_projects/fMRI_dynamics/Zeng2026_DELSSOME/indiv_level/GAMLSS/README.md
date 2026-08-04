# Lifespan E/I-ratio GAMLSS pipeline

Fits a *Generalized Additive Model for Location, Scale and Shape*
([GAMLSS](https://www.gamlss.com/)) to per-participant mean cortical E/I ratios, producing a
normative trajectory over age.

This is the last step of the individual-level pipeline. Get the per-participant E/I ratios first
(see [../README.md](../README.md) §4), then use
`../scripts/CBIG_DELSSOME_indiv_collect_EI.py` to turn them into the input table this pipeline
wants:

```bash
cd .. && export PYTHONPATH=$PWD
python -m scripts.CBIG_DELSSOME_indiv_collect_EI \
    --in-dir results --demogr my_data/demogr.csv --out gamlss_input.csv
```

## Setup

`Rscript` plus the `gamlss` package. The bundled conda environment
([../environment.yml](../environment.yml)) installs both:

```bash
conda env create -f ../environment.yml
conda activate CBIG_DELSSOME_indiv
```

Otherwise, in R: `install.packages(c("gamlss", "gamlss.dist", "gamlss.data", "jsonlite"))`.

## Layout

```
GAMLSS/
├── CBIG_fit_gamlss_model.R          <- main entry point; JSON config + CLI overrides
├── CBIG_apply_gamlss_model.R        <- apply a fitted model to a new cohort
├── CBIG_gamlss_utils.R              <- shared helpers (sourced by fit + apply)
├── CBIG_compare_gamlss_models.R     <- BIC sweep across fitted models
├── CBIG_compare_gamlss_dists.R      <- distribution-only model comparison
├── configs/                         <- JSON configs (distribution x spline DOF)
└── CBIG_generate_fake_gamlss_data.py <- synthesise a fake input table for testing
```

## Input table

`CBIG_fit_gamlss_model.R` reads a CSV with **a header row** and at least these columns:

| Column | Description |
| --- | --- |
| `age_in_years` | Numeric age — fed to the B-spline / polynomial in `μ` and `σ`. |
| `numeric_sex` | `0` (male) or `1` (female) — fixed effect. |
| `site_scanner` | Site × scanner label — modelled as a **random** effect, which doubles as harmonisation. |
| `ei_ratio_mean` | The response variable. |

`CBIG_DELSSOME_indiv_collect_EI.py` produces exactly this, plus `ei_ratio_roi_0 … ei_ratio_roi_67`
so a per-region model can be fitted with `--target`. `CBIG_generate_fake_gamlss_data.py` produces a
synthetic table in the same format.

A single-site cohort still needs the `site_scanner` column; with one level the random effect
degenerates harmlessly and performs no harmonisation.

## Fitting a model

```bash
Rscript CBIG_fit_gamlss_model.R configs/SHASHo2_bs_df4.json \
    --input_data gamlss_input.csv \
    --output_dir outputs/SHASHo2_bs_df4
```

`configs/SHASHo2_bs_df4.json` is the paper's specification: the SHASHo2 (sinh–arcsinh, original
type 2) distribution, cubic B-spline with df=4 on `μ` and df=3 on `σ`, sex as a fixed effect,
site/scanner as a random effect.

Both scripts locate their own helpers, so they can be run from any working directory and with
absolute paths.

The output directory will contain:

```
fitted_model.rds                <- saved GAMLSS object + config
gamlss_outputs.csv              <- input rows + predictions + harmonised values
male_population_centiles.csv    <- normative trajectories for males
female_population_centiles.csv  <- normative trajectories for females
model_summary.txt               <- summary(m) dump
fitting_log.txt                 <- timestamped log of the run
gamlss_console_output.txt       <- RS() iteration log
```

All CLI overrides:

```
Rscript CBIG_fit_gamlss_model.R <config_json> \
    [--input_data <file>] [--target <var>] [--output_dir <dir>] [--distribution <dist>] \
    [--apply_input <file>] [--apply_output <file>] [--seed <int>] [--skip-fit]
```

For example, to model one region instead of the cortical mean:

```bash
Rscript CBIG_fit_gamlss_model.R configs/SHASHo2_bs_df4.json \
    --input_data gamlss_input.csv \
    --target ei_ratio_roi_27 \
    --output_dir outputs/SHASHo2_bs_df4_roi27
```

The configs describe the **model** only — distribution, spline degrees of freedom, output
sub-directory, centiles. The input table is a per-run choice and comes from `--input_data`.

## Applying a fitted model to a new cohort

```bash
Rscript CBIG_apply_gamlss_model.R \
    outputs/SHASHo2_bs_df4/fitted_model.rds \
    new_cohort.csv \
    outputs_apply/new_cohort_harmonised.csv
```

`new_cohort.csv` needs the same columns as the training input. Sites the model has not seen are
mapped to the population-average prediction (no site effect), so the harmonised values stay
well-defined. An `application_log.txt` is written next to the output.

## Model comparison

Once several configs have been fitted, `CBIG_compare_gamlss_models.R` walks every
`fitted_model.rds` under a base directory and prints them ranked by BIC:

```bash
Rscript CBIG_compare_gamlss_models.R <base_dir>
```

`CBIG_compare_gamlss_dists.R` does the same for a `<base_dir>/<distribution>/<target>/` layout and
additionally writes `distribution_comparison.csv` into `<base_dir>`. Both scripts have
`target_var <- "ei_ratio_mean"` as a constant near the top; edit it if you are comparing per-region
models.

## Changing the model

The bundled configs are starting points, not a fixed menu. `configs/` ships `SHASHo2_bs_df4` (the
paper's choice), `jsu_bs_df3/4/5` (Johnson's SU), `gg_bs_df4` (generalised gamma) and
`normal_bs_df4` — but every field in a config is yours to edit.

**A different distribution** takes one line. Copy any config, change `distribution` to the
abbreviation the `gamlss.dist` package uses, and fit:

```json
"distribution": "BCTo",
```

Any distribution documented in
[the gamlss.dist reference manual](https://cran.r-project.org/web/packages/gamlss.dist/gamlss.dist.pdf)
works, provided your response variable lies in its support. Distributions with three or four
parameters also use `age_formula_nu` and `age_formula_tau`; two-parameter families ignore them.

**A different age model** takes one line too. `age_formula_mu` and `age_formula_sigma` (and `_nu` /
`_tau`) are R formula fragments passed straight through, so you can swap the B-spline for a
polynomial, change its degrees of freedom, or hold a parameter constant:

```json
"age_formula_mu":    "bs(age_in_years, df=6, degree=3)",
"age_formula_sigma": "poly(age_in_years, 2)",
"age_formula_nu":    "1",
```

Setting a formula to `"1"` fits that parameter as a constant. Otherwise
`CBIG_fit_gamlss_model.R` appends `+ numeric_sex + random(site_scanner)` to whatever you write, so
the sex fixed effect and the site/scanner random effect are always present; change that in the
script if you need a different design. Fit several configs and rank them with
`CBIG_compare_gamlss_models.R` above.

## Demo

```bash
bash ../../examples/indiv_level/CBIG_DELSSOME_indiv_gamlss_example.sh
```

About 15 seconds. Fits a bundled 600-subject synthetic table and then applies the fitted model to a
second synthetic cohort, exercising both scripts above. Skips cleanly with a `[skip]` notice if
`Rscript` is not on `PATH`.

GAMLSS needs a few hundred participants at minimum before a normative trajectory means anything;
the paper used ~12,000.
