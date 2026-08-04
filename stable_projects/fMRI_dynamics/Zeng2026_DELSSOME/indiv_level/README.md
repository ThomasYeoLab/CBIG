# Individual-level DELSSOME

Estimate a **per-participant, per-region excitation/inhibition (E/I) ratio** from that
participant's own resting-state fMRI, and optionally fit a normative lifespan trajectory across a
cohort.

> The release ships the paper's cost predictor at
> `pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt`. §4 applies it to your data, §5 explains
> what it covers, and §7 — training one of your own — is only for data it does not cover.

Paths in prose are relative to the repository root; commands are run from `indiv_level/`.

This is the individual-level half of the release. The group-level pipeline, which fits
*group-averaged* FC and FCD CDF and also covers the pMFM and Hopf models, is in
[../group_level/](../group_level/).

**Scope**: only the pFIC model (Deco et al., 2014) and only **DK68** — the 68 cortical regions of
the Desikan–Killiany atlas. Hopf and MFM-2013 work with the same architecture (paper Section 2.5)
but are omitted here for clarity.

---

## 1. What you get

| | |
| --- | --- |
| **Input** | one participant's FC matrix + FCD cumulative distribution function |
| **Output** | `EI_ratio.csv` — E/I ratio per cortical region, plus the fitted circuit parameters and the fit quality |
| **Across a cohort** | a GAMLSS normative E/I trajectory over age, with site/scanner harmonisation |
| **Cost** | minutes to hours per participant on one CPU core (§5.2) |

---

## 2. Install

```bash
cd path/to/DELSSOME/indiv_level
conda env create -f environment.yml
conda activate CBIG_DELSSOME_indiv
export PYTHONPATH=$PWD
```

Or `pip install -r requirements.txt`, plus R with the `gamlss` package if you want §6. There is no
`setup.py` — the scripts are run as modules from this directory, which is why `PYTHONPATH` is set
above.

---

## 3. Try the bundled demo first

```bash
bash ../examples/indiv_level/CBIG_DELSSOME_indiv_EI_example.sh
```

A few minutes, no configuration, writes `EI_ratio.csv` for two synthetic subjects — see
[../examples/indiv_level/README.md](../examples/indiv_level/README.md). The data is synthetic, so
the numbers are not meaningful; the point is to confirm your environment works and to see the
input/output shapes.

---

## 4. Run it on your own data

### 4.1 What you need

All matrices are **header-less CSV** — no header row, no index column. N = number of ROIs (68 for
DK68).

| Input | Shape | Notes | Required |
| --- | --- | --- | --- |
| FC | N × N | Pearson correlation between region time series. Symmetric, diagonal 1. | **yes** |
| FCD CDF | bins × 1 | Cumulative distribution of the FCD matrix's off-diagonal entries, over `bins` bins spanning (−1, 1). Non-decreasing, last value 1. Must match `fcd_hist_bins` in the config (10000 by default). | **yes** |
| TR | scalar | Repetition time of your fMRI, **in seconds** | **yes** |
| scan length | scalar | Duration of your scan, **in seconds** (frames × TR) | **yes** |
| SC | N × N | Structural connectivity, zero diagonal | no — defaults to the bundled HCP-YA group SC |
| myelin | N × 1 | Regional myelin (T1w/T2w) profile | no — defaults to the bundled HCP-YA map |
| RSFC gradient | N × 1 | Regional resting-state FC gradient | no — defaults to the bundled HCP-YA map |

`TR` and `scan length` have **no defaults on purpose**. They set the length and sampling of the
simulated BOLD signal, so a wrong value does not raise an error — it quietly produces an E/I ratio
fitted against a mis-specified simulation. Read them off your own acquisition.

SC, myelin and the RSFC gradient default to the HCP-YA group maps because that is what the paper
does — per-subject SC was not available across all 14 datasets. Supply your own with `--sc`,
`--myelin` and `--gradient` if you have them (`--myelin` and `--gradient` must be given together:
they jointly parameterise `wEE`, `wEI` and `sigma`). A non-DK68 parcellation therefore *requires*
your own myelin and gradient maps.

### 4.2 One participant, one command

```bash
python -m scripts.CBIG_DELSSOME_indiv_estimate_EI \
    --fc            my_data/sub-01_FC.csv \
    --fcd-cdf       my_data/sub-01_FCD_CDF.csv \
    --tr            2.0 \
    --scan-length   600 \
    --delssome-ckpt ../pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt \
    --out-dir       results/sub-01
```

`--delssome-ckpt` is the released predictor — it makes the parameter search cheaper and is the
recommended default for DK68 data from healthy participants (§5). Drop the flag to score every
candidate by Euler integration instead: the outputs are identical in format, and the search is
slower. If the checkpoint is not present in your copy of the repository, see
[Data](../README.md#data).

Add `--seed 1` for reproducible output, `--num-epochs` to trade time for search quality (default
100), and `--sc / --myelin / --gradient` to override the bundled maps. `--help` lists everything.

### 4.3 What you get back

| File | Columns |
| --- | --- |
| `EI_ratio.csv` | `region_index, EI_ratio, S_E, S_I, r_E` — the headline result, one row per region |
| `best_parameters.csv` | `region_index, wEE, wEI, sigma, G` — the fitted circuit parameters. `G` is a single global coupling, repeated down the column so the file stays rectangular. |
| `fit_quality.csv` | `stage, corr_loss, l1_loss, ks_loss, total_loss` for the train / val / test stages |
| `run_config.ini` | the fully resolved configuration, for provenance |
| `.work/` | intermediate `.pth` files; safe to delete |

`fit_quality.csv` is worth reading before you trust a result. `corr_loss` is `1 − r` between
simulated and empirical FC, `l1_loss` is their mean absolute difference, and `ks_loss` is the
Kolmogorov–Smirnov distance between simulated and empirical FCD distributions. If the test-stage
`total_loss` is far worse than the train-stage value, the search overfitted and you should increase
`--num-epochs`.

### 4.4 A whole cohort

There is no built-in batch mode; loop in the shell and collect afterwards. Each participant is
independent, so on a cluster submit one job per participant — the command below runs unchanged
inside one.

```bash
for sub_id in $(cat my_subject_list.txt); do
    python -m scripts.CBIG_DELSSOME_indiv_estimate_EI \
        --fc      my_data/${sub_id}_FC.csv \
        --fcd-cdf my_data/${sub_id}_FCD_CDF.csv \
        --tr 2.0 --scan-length 600 \
        --delssome-ckpt ../pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt \
        --subject-id ${sub_id} \
        --out-dir results/${sub_id}
done

# one row per participant, in the format the GAMLSS step expects
python -m scripts.CBIG_DELSSOME_indiv_collect_EI \
    --in-dir results --demogr my_data/demogr.csv --out gamlss_input.csv
```

`--demogr` wants a CSV with a `subject_id` column plus what the GAMLSS formula uses:

```
subject_id,age_in_years,sex,site,scanner
sub-01,34.0,female,WashU,Siemens
```

---

## 5. The released predictor

### 5.1 When it applies

`pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt` is the paper's individual-level predictor.
It was trained on individual participants pooled from **12 healthy-participant datasets**, and it
transfers to datasets that were not in that pool — so a new study does not imply a new predictor,
and §7 is only for the cases below that it does not cover.

| Your data | The released predictor |
| --- | --- |
| DK68, 10000 FCD-CDF bins, healthy participants aged roughly 4–98 | applies — including sites, scanners and datasets it has never seen |
| healthy, but outside roughly 4–98 years | outside the age span it was trained on; untested |
| a clinical or otherwise non-healthy cohort | the training pool was healthy participants only; untested |
| a different parcellation or FCD binning | does not apply. The architecture is tied to both, and `estimate_EI` rejects a mismatched checkpoint rather than failing obscurely later. |

Being outside that envelope degrades the search, not the result (§5.2). To check how well the
predictor holds up on your own data, run a few participants with and without `--delssome-ckpt` and
compare `fit_quality.csv` (§4.3).

### 5.2 What it changes, and how long a run takes

CMA-ES runs in three stages per participant:

| Stage | What it does | How each candidate is scored |
| --- | --- | --- |
| `train` | Searches the parameter space — `--num-epochs` epochs × 100 candidates | Euler integration, **or** the predictor with `--delssome-ckpt` |
| `val` | Re-scores the per-epoch winners at the validation step size | always Euler |
| `test` | Re-simulates the single best parameter set at a finer step size, and reads off the E/I ratio | always Euler |

`--delssome-ckpt` does not change the CMA-ES search itself. It changes how each candidate parameter
set is **evaluated**: instead of Euler-integrating the circuit model to obtain simulated FC and FCD
and then scoring them, the predictor returns the FC+FCD cost in one forward pass. Only the `train`
stage uses it — validation and test still evaluate by Euler simulation (paper Sec. 4.2), so the
E/I ratio you read off is never a neural-network estimate, and a predictor that transfers poorly
costs you search efficiency rather than correctness.

DELSSOME therefore does not make this interactive. Budget minutes to hours per participant either
way, dominated by the test stage: its Euler step size `dt_test` is 10× finer than the `dt_train`
used during the search, so each test simulation costs 10× more steps.

---

## 6. Normative lifespan trajectories (GAMLSS)

Once you have a `gamlss_input.csv` from §4.4:

```bash
Rscript GAMLSS/CBIG_fit_gamlss_model.R GAMLSS/configs/SHASHo2_bs_df4.json \
    --input_data gamlss_input.csv --output_dir results/gamlss
```

`SHASHo2_bs_df4.json` is the paper's specification: the SHASHo2 (sinh–arcsinh, original type 2)
distribution, cubic B-spline with df=4 on µ and df=3 on σ, sex as a fixed effect, and site/scanner
as a random effect — which doubles as site harmonisation, so no separate ComBat step is used. See
[GAMLSS/README.md](GAMLSS/README.md) for the table format, how to switch distribution or change the
µ/σ formulas, model comparison by BIC, and applying a fitted model to a new cohort.

GAMLSS needs a few hundred participants at minimum; the paper used ~12,000.

---

## 7. Advanced: training your own predictor

**Read §5.1 first — most cohorts do not need this section.** Train your own predictor only when
your parcellation, FCD binning or cohort falls outside what the released one covers — and even
then, only for a large cohort. Below roughly 20 participants, training a predictor costs more than
the evaluations it would save, so just run pure Euler (drop `--delssome-ckpt`).

```bash
# 1. Generate training data: Euler CMA-ES on a pilot subset (the expensive part).
#    Each subject's 100 epochs x 100 candidates give 10,000 (parameters, cost) pairs.
for sub_id in $(cat pilot_subjects.txt); do
    python -m scripts.CBIG_DELSSOME_indiv_cmaes \
        --ds-name MYSTUDY --sub-id ${sub_id} --engine euler \
        --mode all --num-epochs 100 --trial 1 --seed-idx 1
done

# 2. Train the predictor (Optuna hyper-parameter sweep; use a GPU).
python -m scripts.CBIG_DELSSOME_indiv_predictor_trainer \
    --trial 1 --seed 1 --ds-names MYSTUDY \
    --max-epochs 15 --n-trials 50 --walltime 120:00:00 \
    --batch-size 2048 --n-sample-per-phase train40_val20_test20

# 3. (optional) Evaluate it on the held-out split.
python -m scripts.CBIG_DELSSOME_indiv_predictor_tester --ckpt-path <path-to.ckpt>

# 4. Apply it to the rest of the cohort -- via estimate_EI (§4.2) or:
python -m scripts.CBIG_DELSSOME_indiv_cmaes \
    --ds-name MYSTUDY --sub-id ${sub_id} --engine delssome \
    --mode all --num-epochs 100 --trial 2 --seed-idx 1 \
    --delssome-ckpt <path-to.ckpt>
```

For reference, the paper's individual-level predictor took ~65 h on one RTX 3090 for 50 Optuna
trials, and settled on a transformer with hidden dimension 128, 16 attention heads and 5 layers.

### 7.1 The directory tree these scripts expect

Steps 1, 2 and 4 use `CBIG_DELSSOME_indiv_cmaes.py` and the trainer, which — unlike `estimate_EI` —
organise everything through a directory tree rather than explicit file paths, because predictor
training has to walk many subjects at once.

Four environment variables name the corners of that tree. Export them before running anything in
this section; each one falls back to a location inside `indiv_level/`, which is only convenient for
the bundled demo:

| Variable | What lives there | Default |
| --- | --- | --- |
| `DELSSOME_PROJECT_DIR` | the `indiv_level/` tree itself; the other three defaults hang off it | the `indiv_level/` folder containing this README |
| `DELSSOME_DATA_DIR` | your inputs: per-subject FC and FCD CDF, demographics, the shared group maps | `$DELSSOME_PROJECT_DIR/examples/data` |
| `DELSSOME_LOG_DIR` | outputs: split lists and per-subject CMA-ES results | `$DELSSOME_PROJECT_DIR/examples/logs` |
| `DELSSOME_CONFIG_DIR` | the `configs/` tree of `.ini` files | `$DELSSOME_PROJECT_DIR/configs` |

The layout they expect:

```
$DELSSOME_DATA_DIR/<ds_name>/
    ├── FC/<atlas>/<sub_id>_bld000.csv        # one FC matrix per subject
    ├── FCD_cdf/<atlas>/<sub_id>_bld000.csv   # one FCD CDF per subject
    └── demogr/demogr.csv                     # subject_id, age_in_years, sex, site, scanner, TR, scan_length
$DELSSOME_DATA_DIR/HCP-YA/pFIC_input/<atlas>/
    ├── SC_train_group_dl_ds.csv              # group SC, shared by every dataset
    ├── myelin.csv
    └── rsfc_gradient.csv

$DELSSOME_LOG_DIR/<ds_name>/<target>/<stage>/sub_list.txt                 # who is in each split
$DELSSOME_LOG_DIR/<ds_name>/<target>/<stage>/trial<T>/seed<S>/<sub_id>/   # results
```

`bld000` in those filenames is the across-run average. If a subject has several fMRI runs you may
instead supply them separately as `<sub_id>_bld001.csv` … `<sub_id>_bld006.csv`, and the loader
Fisher-averages the FC matrices (and plain-averages the FCD CDFs) for you.

The `HCP-YA/pFIC_input/` directory is required **even when your dataset is not HCP-YA** — it is
where the shared group SC / myelin / gradient maps live.

Each dataset also needs its own config at
`$DELSSOME_CONFIG_DIR/model/pfic/<atlas>/indiv/config_<ds_name>.ini`; copy
[`config_template.ini`](configs/model/pfic/DK68/indiv/config_template.ini) and set your scan timing.

---

## 8. Glossary

Terms that appear in the code and in the paper, and are easy to misread:

- **stage** — one of the three phases CMA-ES runs for a single participant: `train`, `val`, `test`,
  each with its own Euler step size `dt_train` / `dt_val` / `dt_test` (§5.2). `estimate_EI` runs all
  three for you, so you never choose one.
- **split** — how *participants* are divided into train / validation / test for predictor training
  (§7). In the paper's experiment, HCP-YA's 1029 participants split 680 / 180 / 169.
- So **`train` / `val` / `test` mean two different things**: three stages of one participant's
  CMA-ES run, or three splits of a cohort. Which one is meant depends on what is being divided.
- **group** — one bootstrap sample of 50 participants drawn within a split, from which
  group-averaged SC / FC / FCD-CDF are computed. The paper uses 64 / 14 / 13 groups for
  train / val / test. Only relevant to the group-level pipeline; see
  [../group_level/README.md](../group_level/README.md).
- **trial** / **seed_idx** / **seed** — `trial` and `seed_idx` are folder labels only, used to keep
  multiple runs apart; `seed` is the actual RNG seed. `estimate_EI` exposes only `--seed`.
- **target** — the log-folder name identifying a kind of run: `indiv_Euler-pfic`,
  `train_indiv_DELSSOME-pfic`, `apply_indiv_DELSSOME-pfic`.
- **FC+FCD cost** — the quantity CMA-ES minimises and DELSSOME predicts:
  `corr_loss + l1_loss + ks_loss` (§4.3).

---

## 9. Code layout

```
indiv_level/
├── scripts/
│   ├── CBIG_DELSSOME_indiv_estimate_EI.py       <- START HERE: your CSVs -> EI_ratio.csv
│   ├── CBIG_DELSSOME_indiv_collect_EI.py        <- per-subject results -> GAMLSS table
│   ├── CBIG_DELSSOME_indiv_cmaes.py             <- lower-level CMA-ES driver (--engine euler|delssome)
│   ├── CBIG_DELSSOME_indiv_predictor_trainer.py <- train a cost predictor (Optuna)
│   └── CBIG_DELSSOME_indiv_predictor_tester.py  <- evaluate a trained predictor
├── DELSSOME_indiv/
│   ├── CBIG_constants.py        <- paths + dataset / atlas constants
│   ├── CBIG_file_utils.py       <- path builders, config loading, empirical-data loaders
│   ├── CBIG_FC_FCD_utils.py     <- FC, FCD, and the KS / L1 / corr cost components
│   ├── CBIG_misc_utils.py       <- correlation / parameterisation helpers, torch defaults
│   ├── CBIG_dl_train_utils.py   <- Optuna + Lightning plumbing
│   ├── models/
│   │   ├── CBIG_dynamic_model.py       <- pFIC: Euler integration + Balloon-Windkessel
│   │   └── CBIG_DELSSOME_predictor.py  <- the transformer cost predictor (paper Fig. 9)
│   ├── datasets/CBIG_dl_dataset.py     <- DataModules over per-subject Euler dumps
│   └── optimizers/CBIG_cmaes.py        <- CmaesTrainer / ModelValidator / ModelTester
├── configs/model/pfic/DK68/indiv/
│   ├── config_HCP-YA.ini        <- released defaults
│   └── config_template.ini      <- copy this for your own dataset
└── GAMLSS/                      <- normative lifespan fit (R)
```

Implementation notes:

- **No hard-coded paths.** Everything derives from the `DELSSOME_*` environment variables, which
  default to locations under `indiv_level/`.
- **Extending to other circuit models**: register a class in
  `DELSSOME_indiv.models.CBIG_dynamic_model.AVAIL_MODELS`.
- **No ComBat harmonisation**: this release works directly with raw per-subject FC/FCD;
  harmonisation happens inside GAMLSS via the `random(site_scanner)` term (paper Sec. 2.7).
- **No firing-rate gate and no within-range classifier**: following the paper (Sec. 4.2), the
  FC+FCD cost is the sole selection criterion.

