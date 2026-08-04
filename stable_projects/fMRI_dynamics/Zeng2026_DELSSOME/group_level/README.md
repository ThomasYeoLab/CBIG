# DELSSOME — group-level pipeline

> **Looking for per-participant E/I ratios?** They come from the individual-level pipeline in
> **[../indiv_level/](../indiv_level/)**, not from this folder. That is what most people want:
> `indiv_level/scripts/CBIG_DELSSOME_indiv_estimate_EI.py` takes one participant's FC and FCD CDF
> and writes a per-region E/I ratio CSV. Start at
> [../examples/indiv_level/](../examples/indiv_level/).

This folder estimates **one E/I ratio map for a group** of participants, by fitting group-averaged
FC and FCD CDF with CMA-ES. It is also the pipeline that trains the DELSSOME cost predictor, and the
only one that covers the pMFM and Hopf models.

Instead of Euler-integrating the circuit model for every candidate
parameter set, a pretrained deep-learning predictor scores each candidate in a single forward pass —
the bundled pFIC demo finishes in well under a minute where the equivalent Euler run takes long,
and the gap widens with the number of epochs. Validation and testing still evaluate by Euler
simulation, so the parameter set you take away has been verified by real simulation.

> The release ships one per circuit model under
> `pretrained_models/group_level/`, and every bundled config already points at the right one, so §2
> works out of the box. §3 describes them; §4 — training your own — is for a different
> parcellation or out-of-distribution data.

Replication-only utilities (SBI benchmarks, deep-learning dataset generators) and full-dataset
config files are kept under [../replication/](../replication/) and are not part of this release.

---

## 1. Run the bundled demos first

```bash
cd path/to/DELSSOME/group_level
conda env create -f ../replication/config/CBIG_DELSSOME_python_env.yml
conda activate CBIG_DELSSOME

cd ../examples/group_level/scripts
bash CBIG_DELSSOME_CMAES_pFIC_DL_example.sh   # pFIC in DELSSOME mode, < 1 min
bash CBIG_DELSSOME_CMAES_pFIC_example.sh      # the same fit in pure Euler mode, ~2 min
```

Both write the same kind of output; the first one is the fast path.

Four demos in total, with run times and expected outputs, are described in
[../examples/group_level/README.md](../examples/group_level/README.md). Their configs in
[../examples/group_level/config/](../examples/group_level/config/) are the easiest starting point
for your own.

---

## 2. Invert a circuit model

### 2.1 The command

All scripts here import sibling modules by bare name, so **run them with this directory as the
working directory** — the example shell scripts do this for you.

```bash
cd path/to/DELSSOME/group_level
python -u CBIG_DELSSOME_CMAES_main_pFIC.py \
    --mode        train_val_EI \
    --config_path ../examples/group_level/config/example_CMAES_pFIC_DELSSOME.ini \
    --input_dir   my_data/ \
    --output_dir  results/ \
    --n_seed 1 --num_epochs 100
```

All three scripts run in **`--train_mode DELSSOME` by default**: each candidate parameter set is
scored by the cost predictor named by `predictor_save_path` in the config — the released one, unless
you change it (§3). Pass `--train_mode Euler` to integrate the model directly instead; that needs no
predictor, and it is what the validation and test phases use in either case.

`CBIG_DELSSOME_CMAES_main_pMFM.py` and `CBIG_DELSSOME_CMAES_main_Hopf.py` take the same interface,
minus the `EI` modes, and each has a ready config —
[`example_CMAES_pMFM.ini`](../examples/group_level/config/example_CMAES_pMFM.ini) and
[`example_CMAES_Hopf.ini`](../examples/group_level/config/example_CMAES_Hopf.ini). One thing to set
explicitly for a real run: `--num_epochs` defaults to `2` in both scripts (100 in pFIC). pMFM is the
one model without a bundled demo script; run it with the command above, swapping in the pMFM entry
point and config.

### 2.2 Arguments

| Argument | Meaning |
| --- | --- |
| `--config_path` | path to the `.ini` config (see [../examples/group_level/config/](../examples/group_level/config/)) |
| `--input_dir` | directory holding the input files (§2.3) |
| `--output_dir` | directory for the outputs (§2.4) |
| `--train_mode` | `DELSSOME` (cost predictor, the default) or `Euler` (direct simulation) |
| `--mode` | which phases to run — see below |
| `--num_epochs` | total CMA-ES training epochs (default 100 for pFIC, 2 for pMFM and Hopf) |
| `--n_seed` | number of random initializations, run **sequentially** in one process (default 1). See §2.5 before raising it. |
| `--start_seed_nbr` | which seed number this process runs (default 1); names the `seed<n>/` output folder |
| `--seed` | fix the RNG seed(s) for replication. If unset, a seed is chosen and saved for you. |
| `--trials` | how many trials are used to select new seeds (default 10; forced to 1 when `--seed` is given) |
| `--next_epoch` | resume CMA-ES from this epoch index (0-based) |
| `--GPU_index` | accepted but not wired up: the CMA-ES entry points run on CPU. Predictor training and fine-tuning (§4) use the GPU when one is available. |
| `--val_by_dl` | pFIC only: validate with the predictor as well. Off by default — the paper validates by Euler simulation; see §3. |

`--mode` selects the phases:

| `--mode` | Runs |
| --- | --- |
| `train` | the CMA-ES search that generates candidate parameter sets |
| `validation` | re-score existing candidates on the validation set |
| `test` | evaluate the best validation parameter set on the test set |
| `EI` | estimate the E/I ratio from the parameters with the lowest validation loss (pFIC only) |
| `train_val_test` | training, validation and testing in sequence |
| `train_val_EI` | training, validation and E/I estimation in sequence (pFIC only) |

### 2.3 Input data format

`--input_dir` must contain, for each split (`train`, `validation`, `test`):

* `SC_<split>.csv`: structural connectivity matrix (N × N), N = number of ROIs
* `FC_<split>.csv`: functional connectivity matrix (N × N)
* `FCD_CDF_<split>.mat`: FCD CDF matrix (MATLAB format, field `FCD_CDF`)
* `myelin.csv`: myelin map (N × 1)
* `rsfc_gradient.csv`: RSFC gradient (N × 1)

For E/I estimation, additionally `SC_EI.csv`, `FC_EI.csv` and `FCD_CDF_EI.mat`.

### 2.4 Output file format

**Training** (`<output_dir>/train/seed<s>/`)

* `param_save_epoch<e>.pth`:
    * `parameter`: parameter sets ((3N+1) × M for pFIC), N = number of ROIs, M = total sampled
      parameter sets. Arranged as wEE (N), wEI (N), G (1, long-range coupling), sigma (N, neuronal
      noise).
    * `valid_param_list`: indices (K × 1) of parameters whose excitatory firing rates are within
      the allowed range (e.g. 2.7–3.3 Hz).
    * `FC_FCD_loss`: costs (K × 3) — FC correlation loss, FC L1 loss, FCD KS loss.
* `final_state_pFIC.pth` (pFIC) / `final_state.pth` (pMFM, Hopf): `seed` for replication,
  `random_state` for resuming, and the CMA-ES state `m`, `sigma`, `cov`, `p_sigma`, `p_c`.

**Validation** (`<output_dir>/validation/seed<s>/`)

* `best_params_all_epochs.pth`: per-epoch best parameter set with `corr_loss`, `l1_loss`, `ks_loss`.

**Test** (`<output_dir>/test/seed<s>/`)

* `test_results.pth`: `parameter` (the set with the lowest validation loss across all epochs),
  `val_loss`, and the test-set `corr_loss`, `l1_loss`, `ks_loss`.

**E/I ratio** (`<output_dir>/EI_ratio/`)

* `seed<s>.pth`: `ei_ratio` (N × 1), `s_e_ave` and `s_i_ave` (simulated excitatory and inhibitory
  activity), `parameter`, `seed`, `time_series`.

### 2.5 Running several seeds

CMA-ES is stochastic, so a real analysis usually repeats the optimization from several random
initializations and keeps the best. Do **not** recommend doing that by raising `--n_seed`: a single process
works through the seeds sequentially, so `--n_seed 5` simply takes five times as long.

The seeds are completely independent, so run one process per seed and let the operating system
place them on separate CPU cores:

```bash
for seed_nbr in 1 2 3 4 5; do
    python -u CBIG_DELSSOME_CMAES_main_pFIC.py \
        --mode train_val_test --config_path <config_path> \
        --input_dir <input_dir> --output_dir <output_dir> \
        --n_seed 1 --start_seed_nbr ${seed_nbr} &
done
wait
```

`--start_seed_nbr` selects which seed a process runs, and its results land under `seed<seed_nbr>/`
in the shared `--output_dir`, so the parallel runs do not collide. On a cluster, submit the same
command as one job per seed instead of backgrounding it. The same applies to pMFM and Hopf.

---

## 3. The released cost predictors

One predictor per circuit model, all used automatically in the default `DELSSOME` mode:

| File | Used by | Bundled config |
| --- | --- | --- |
| `pretrained_models/group_level/pFIC_DK68_predictor.pth` | `CBIG_DELSSOME_CMAES_main_pFIC.py` | `example_CMAES_pFIC_DELSSOME.ini`, `example_CMAES_pFIC.ini` |
| `pretrained_models/group_level/pMFM_DK68_predictor.pth` | `CBIG_DELSSOME_CMAES_main_pMFM.py` | `example_CMAES_pMFM.ini` |
| `pretrained_models/group_level/Hopf_DK68_predictor.pth` | `CBIG_DELSSOME_CMAES_main_Hopf.py` | `example_CMAES_Hopf.ini` |

Every bundled config already points at the right one, in its `[Deep learning Model Constants]`
section — for example
[`example_CMAES_pFIC_DELSSOME.ini`](../examples/group_level/config/example_CMAES_pFIC_DELSSOME.ini):

```ini
[Deep learning Model Constants]
predictor_save_path = ../pretrained_models/group_level/pFIC_DK68_predictor.pth
d_transformer = 64
```

The path is resolved relative to the working directory, which is this folder (§2.1). `d_transformer`
must match the width the checkpoint was trained with — 64 for all three released predictors. If the
file is not present in your copy of the repository, see [Data](../README.md#data).

**What they cover.** All three are **DK68** (68 cortical regions) and were trained on group-averaged
connectivity bootstrapped from HCP-YA (§4.2), so any DK68 group-averaged data runs against them
directly. The paper applies them across cohorts; §4 is there for the one case they cannot cover, a
different parcellation, because the architecture is sized by `n_regions`.

**Your results are always simulation-verified.** By design, only the training phase uses the
predictor — validation and testing evaluate by Euler simulation, as recommended in the paper — so
the parameter set you end up with has been scored by the circuit model itself. The predictor buys
search speed; correctness is still established by simulation. That is what makes DELSSOME safe to
use on a new cohort: at worst the search is less efficient, never silently wrong. (`--val_by_dl`
moves the validation phase onto the predictor too; the paper does not use it.)

---

## 4. Advanced: training your own predictor

This is the path for training a new predictor, and it is also how the predictors in §3 were made, so
the same recipe reproduces the paper's training pipeline end to end.

### 4.1 The scripts

* `CBIG_DELSSOME_predictor_trainer.py` — trains the FC+FCD cost predictor from scratch. One
  invocation always trains and then tests.
    * Arguments: `--config_path`, `--seed`
    * Config sections: `[data preparation]`, `[predictor training]`, `[predictor model]`
    * Demo: `CBIG_DELSSOME_predictor_train_example.sh`, about a minute on a tiny subset.

* `CBIG_DELSSOME_finetuner.py` — fine-tunes or retrains a predictor. One script covers pFIC, pMFM
  and Hopf via two `[predictor training]` flags:
    * `finetuning` (bool): load `pre_trained_model_path` before training (pMFM) vs. train from
      scratch (pFIC/Hopf, which the authors found works better).
    * `freeze_transformer` (bool): freeze the transformer layers during fine-tuning (pMFM setup).
    * Arguments: `--config_path`, `--seed`
    * Config: [`example_finetune_pFIC.ini`](../examples/group_level/config/example_finetune_pFIC.ini).
      Run the predictor-training demo above first, so there is a checkpoint to load.

### 4.2 Training data: `splits`, `groups` and `bootstrap_mats`

Two words that are easy to confuse, and that several config keys depend on:

* A **split** partitions *participants* into train / validation / test. The paper divides HCP-YA's
  1029 participants into 680 / 180 / 169. Splits are made first so that group-level statistics
  never leak across them.
* A **group** is one sample of 50 participants drawn *within* a split, from which one
  group-averaged SC, FC and FCD-CDF triplet is computed.

680 participants allow at most 13 non-overlapping groups of 50, too few to capture the variability
of group-level connectivity, so the paper samples 50 participants repeatedly — 64 / 14 / 13 times
in the train / validation / test splits. Each triplet is then inverted with 100 epochs of Euler
CMA-ES at 100 candidates per epoch, yielding 10,000 `(parameters, ground-truth FC+FCD cost)` pairs
per group: 640,000 / 140,000 / 130,000 samples in total.

`data/bootstrap_mats/{train,val,test}.mat` holds those group-averaged matrices, one MATLAB file per
split, with three variables:

| Variable | Shape | Contents |
| --- | --- | --- |
| `sc_groups` | N × N × n_groups | group-averaged structural connectivity |
| `fc_groups` | N × N × n_groups | group-averaged functional connectivity |
| `fcd_cdf_groups` | bins × n_groups | group-averaged FCD cumulative distribution |

To build your own, for each group of participants:

1. **SC** — threshold likely false positives first: if fewer than 50% of participants have a
   non-zero value in an entry, zero that entry in every participant's matrix
   (de Reus & van den Heuvel, 2013). Average the remaining non-zero values across participants,
   log-transform, zero the diagonal, and normalise so the maximum is 0.02.
2. **FC** — Fisher-z average the participants' correlation matrices, then transform back.
3. **FCD CDF** — compute each participant's FCD matrix with a ~60 s sliding window, histogram its
   off-diagonal entries over `bins` bins spanning (−1, 1), take the cumulative sum, normalise so
   the last value is 1, and average across participants.

`data/param_dataset/<split>/group<N>/` then holds the CMA-ES output for that group —
`param_save_epoch{0..99}.pth` plus `final_state_pFIC.pth` — which is what the predictor actually
trains on. `[data preparation]` in the config selects how much of it to load:
`n_training_group` / `n_val_group` / `n_test_group` (groups per split) and `data_epochs` (epochs per
group). The bundled example reduces these from 64 / 14 / 13 / 100 to 2 / 2 / 2 / 10 so it finishes
in about a minute.

---

## 5. Code layout

```
group_level/
├── CBIG_DELSSOME_CMAES_main_{pFIC,pMFM,Hopf}.py  <- START HERE: CMA-ES entry points
├── CBIG_DELSSOME_{pFIC,pMFM,Hopf}_optimizer.py   <- the core optimization logic per model
├── CBIG_DELSSOME_predictor_trainer.py            <- train a cost predictor
├── CBIG_DELSSOME_finetuner.py                    <- fine-tune / retrain one
├── CBIG_DELSSOME_pFIC_get_EI.py                  <- E/I estimation; a library function
│                                                    imported by the pFIC main script, no CLI
├── models/
│   ├── CBIG_DELSSOME_pFIC.py     <- pFIC: Euler integration, neural -> BOLD, cost function
│   ├── CBIG_DELSSOME_pMFM.py     <- parametric Mean-Field Model
│   ├── CBIG_DELSSOME_Hopf.py     <- Hopf bifurcation model
│   └── CBIG_DELSSOME_model.py    <- the DELSSOME cost-predictor architecture
├── process/
│   ├── CBIG_DELSSOME_data_process.py  <- retrieves the inputs for predictor training
│   └── CBIG_DELSSOME_datasets.py      <- PyTorch Datasets for parameter sets, losses, connectivity
└── utils/CBIG_DELSSOME_pFIC_utils.py  <- FC/FCD computation, parameter conversion, correlation
```

## 6. Hyperparameters

See [pFIC_hyperparameter_readme.md](pFIC_hyperparameter_readme.md) for detailed explanations of
every pFIC hyperparameter, grouped by how likely you are to need to change it.
