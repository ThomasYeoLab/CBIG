# Individual-level DELSSOME — examples

Two demos, both on bundled synthetic data:

```bash
cd /path/to/DELSSOME
conda activate CBIG_DELSSOME_indiv

bash examples/indiv_level/CBIG_DELSSOME_indiv_EI_example.sh       # FC + FCD CDF -> E/I ratio
bash examples/indiv_level/CBIG_DELSSOME_indiv_gamlss_example.sh   # cohort of E/I ratios -> lifespan curve
```

> The bundled data is **synthetic** — random matrices with the right shapes, generated from a fixed
> seed (see [../../README.md](../../README.md), "Example data provenance"). The E/I values these
> demos print are therefore meaningless as neuroscience. They demonstrate the mechanics and detect
> regressions. Real numbers need real FC and FCD CDFs.

---

## 1. `CBIG_DELSSOME_indiv_EI_example.sh` — E/I ratio for one subject

This is the whole user-facing workflow. For each of two demo subjects it runs

```bash
python -m scripts.CBIG_DELSSOME_indiv_estimate_EI \
    --fc          input/sub-000_FC.csv \
    --fcd-cdf     input/sub-000_FCD_CDF.csv \
    --sc          input/SC_group_DK68.csv \
    --myelin      input/myelin_DK68.csv \
    --gradient    input/rsfc_gradient_DK68.csv \
    --tr          0.72 \
    --scan-length 60 \
    --config      config/example_indiv_pFIC.ini \
    --num-epochs  3 \
    --seed        123456 \
    --out-dir     output/euler/sub-000
```

then assembles the results into a GAMLSS input table with
`CBIG_DELSSOME_indiv_collect_EI.py`, and — if the pretrained individual-level predictor is
present at `../../pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt` — repeats the first subject
with `--delssome-ckpt` to show the accelerated variant. Without that file the demo still runs;
pure Euler needs no predictor.

**Inputs**: `input/` (four subjects, though the demo uses two by default).
**Outputs**:

| File | Contents |
| --- | --- |
| `output/euler/<sub_id>/EI_ratio.csv` | `region_index, EI_ratio, S_E, S_I, r_E` — one row per cortical region |
| `output/euler/<sub_id>/best_parameters.csv` | `region_index, wEE, wEI, sigma, G` (`G` is global, repeated) |
| `output/euler/<sub_id>/fit_quality.csv` | FC+FCD cost per stage (`corr_loss`, `l1_loss`, `ks_loss`, `total_loss`) |
| `output/euler/<sub_id>/run_config.ini` | the fully resolved config, for provenance |
| `output/euler/gamlss_input_from_demo.csv` | the collected table, in the format the GAMLSS step wants |
| `output/delssome/<sub_id>/` | same files from the DELSSOME-accelerated run, when a predictor is available |

**Check**: `CBIG_DELSSOME_indiv_check_EI_example_results.py`, tolerance 1e-6 on the summed absolute
difference against `reference_output/` — effectively bit-exact, since CMA-ES is seeded.

**Overridable**: `SUB_IDS`, `NUM_EPOCHS`, `SEED`, `TR`, `SCAN_LENGTH`, `DELSSOME_CKPT`, `SKIP_CHECK`.

### Why the demo is fast, and what that costs

`config/example_indiv_pFIC.ini` shrinks the simulation so the demo finishes in minutes rather than
about an hour: 60 s of simulated BOLD instead of 864 s, a 10-frame FCD window instead of 83, a 6 s
burn-in instead of 72, and — the big one — `dt_test = 0.006` instead of `0.0005`, which is 12×
coarser. A real analysis should drop `--config` entirely and use the released defaults in
`indiv_level/configs/model/pfic/DK68/indiv/config_HCP-YA.ini`. This is also why the demo's E/I
values should not be compared with published numbers.

---

## 2. `CBIG_DELSSOME_indiv_gamlss_example.sh` — normative lifespan trajectory

Fits a Generalized Additive Model for Location, Scale and Shape to a bundled **600-subject
synthetic** lifespan E/I table, then applies the fitted model to a second synthetic cohort to
demonstrate the site-harmonisation workflow.

It deliberately does **not** consume the output of demo 1: two subjects cannot support a normative
model. `input/gamlss_input.csv` shows the exact format that
`CBIG_DELSSOME_indiv_collect_EI.py` produces for your own cohort, so the two demos join up as soon
as you have a few hundred subjects.

Default config is `indiv_level/GAMLSS/configs/SHASHo2_bs_df4.json` — the paper's choice: the SHASHo2
(sinh–arcsinh, original type 2) distribution, cubic B-spline with df=4 on µ and df=3 on σ, sex as a
fixed effect, site/scanner as a random effect (which doubles as harmonisation).

**Outputs**: `output/gamlss/<config>/` (fitted model, centiles, summary) and
`output/gamlss/apply_demo/` (harmonised values for the new cohort).
**Check**: `CBIG_DELSSOME_indiv_check_gamlss_example_results.py`, tolerance 1e-4 absolute — looser
than demo 1 because the GAMLSS convergence path depends on the BLAS and the R/gamlss build. Only
the default config ships a reference snapshot; with any other `CONFIG_NAME` the fit still runs and
the comparison is skipped.
**Overridable**: `CONFIG_NAME` (any JSON under `indiv_level/GAMLSS/configs/`), `SKIP_CHECK`.
**Skipped automatically** when `Rscript` is not on `PATH`, with a pointer to
[../../indiv_level/GAMLSS/README.md](../../indiv_level/GAMLSS/README.md).

---

## 3. Bundled input format

This is exactly the format your own data needs; see
[../../indiv_level/README.md](../../indiv_level/README.md) §4 for the full specification.

| File | Shape | Format |
| --- | --- | --- |
| `input/sub-00N_FC.csv` | 68 × 68 | header-less CSV, symmetric, diagonal 1 |
| `input/sub-00N_FCD_CDF.csv` | 10000 × 1 | header-less CSV, non-decreasing, last value 1 |
| `input/SC_group_DK68.csv` | 68 × 68 | header-less CSV, zero diagonal |
| `input/myelin_DK68.csv` | 68 × 1 | header-less CSV |
| `input/rsfc_gradient_DK68.csv` | 68 × 1 | header-less CSV |
| `input/demogr.csv` | 4 rows | **has a header**: `subject_id, age_in_years, sex, site, scanner, TR, scan_length` |
| `input/gamlss_input.csv` | 600 rows | **has a header**: `subject_id, age_in_years, numeric_sex, site_scanner, ei_ratio_mean, ei_ratio_roi_0..67` |

---

## 4. Re-running cleanly

```bash
rm -rf examples/indiv_level/output
bash examples/indiv_level/CBIG_DELSSOME_indiv_EI_example.sh
```

`input/` and `reference_output/` are read-only release artefacts — do not edit them by hand.

The reference comparison exists to catch regressions in the shipped configuration; it is not
something you need to maintain for your own runs. Set `SKIP_CHECK=1` to run a demo without it,
which is what you want once you have started changing the config — the pipeline still writes
everything to `output/`.
