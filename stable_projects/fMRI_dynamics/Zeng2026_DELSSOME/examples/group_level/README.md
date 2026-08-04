# Group-level DELSSOME — examples

Demo input data, configuration files and reference outputs for the group-level pipeline
([../../group_level/](../../group_level/)). For the individual-level demos see
[../indiv_level/](../indiv_level/).

The reference outputs were generated with **CPU only**. To replicate them exactly, use the same
Python package versions — pinned in
[../../replication/config/CBIG_DELSSOME_python_env.yml](../../replication/config/CBIG_DELSSOME_python_env.yml):

```
conda env create -f ../../replication/config/CBIG_DELSSOME_python_env.yml
conda activate CBIG_DELSSOME
```

The installation may take a few minutes.

The pFIC and Hopf examples run out of the box from the inputs committed under `input/`. The
**DELSSOME** (`..._pFIC_DL_...`) and **predictor-training** examples additionally need the
pretrained cost predictor / training matrices, which exceed this repository's 1 MB per-file limit
and are not committed here — obtain them from the stand-alone public repository
(**https://github.com/TianCZeng/DELSSOME**) and place them at `../../pretrained_models/` and
`../../data/`. For CBIG unit tests they are staged automatically from `$CBIG_TESTDATA_DIR`
(see `../../unit_tests/README.md`).

## Usage

The four demo scripts live in `scripts/`. Each one activates the `CBIG_DELSSOME` conda environment,
runs the pipeline, and then invokes its checker to diff the result against `reference_output/` — so
every example doubles as a regression test, and a non-zero exit means a regression.

All run times were measured on a single core of an **AMD EPYC 9555** (64-core, 3.20 GHz, 256 MB
cache, 360 W, 12-channel DDR5-6000); your hardware will differ.

* `CBIG_DELSSOME_CMAES_pFIC_example.sh`

  ```
  cd /path/to/examples/group_level/scripts
  bash CBIG_DELSSOME_CMAES_pFIC_example.sh
  ```

  Runs **pFIC** parameter optimization with CMA-ES in **Euler mode** (no deep-learning cost
  predictor needed) and saves results under `examples/group_level/output/pfic/`, which is created on
  first run. See [../../group_level/README.md](../../group_level/README.md) for the output-file
  format.

  Typical run time: about **2 minutes** on a single CPU.

* `CBIG_DELSSOME_CMAES_pFIC_DL_example.sh`

  ```
  cd /path/to/examples/group_level/scripts
  bash CBIG_DELSSOME_CMAES_pFIC_DL_example.sh
  ```

  Runs **pFIC** parameter optimization with CMA-ES in **DELSSOME mode**, where a pretrained
  deep-learning cost predictor scores each candidate parameter set in one forward pass instead of
  Euler-integrating the model. Requires
  `pretrained_models/group_level/pFIC_DK68_predictor.pth` (DK68 = the 68 cortical regions of the
  Desikan–Killiany atlas). Results go to `examples/group_level/output/dl/`.

  Typical run time: **under a minute** on a single CPU.

* `CBIG_DELSSOME_CMAES_Hopf_example.sh`

  ```
  cd /path/to/examples/group_level/scripts
  bash CBIG_DELSSOME_CMAES_Hopf_example.sh
  ```

  Runs **Hopf bifurcation model** parameter optimization with CMA-ES in **Euler mode** and saves
  results under `examples/group_level/output/hopf/`.

  Typical run time: about **5 minutes** on a single CPU.

* `CBIG_DELSSOME_predictor_train_example.sh`

  ```
  cd /path/to/examples/group_level/scripts
  bash CBIG_DELSSOME_predictor_train_example.sh
  ```

  Trains a DELSSOME cost predictor from scratch on a small subset of the data — **2 groups per
  split** and **10 data epochs per group**, reduced from the paper's 64 / 14 / 13 groups and 100
  epochs. The trained model lands in `examples/group_level/output/predictor/predictor_1/`.

  Here a **split** is the train / validation / test partition of *participants*, and a **group** is
  one sample of 50 participants drawn within a split, from which group-averaged SC, FC and FCD-CDF
  are computed. [../../group_level/README.md](../../group_level/README.md) explains both terms, and
  documents how to prepare `data/bootstrap_mats/` for your own cohort.

  Typical run time: about **1 minute** on a single CPU.

  Unlike the three CMA-ES checkers, which require agreement to 1e-6 on the summed absolute
  difference, this example's checker allows 1e-2 on `mse_mean`: the Optuna and Lightning RNG paths
  are not bit-reproducible across BLAS versions or CPU/CUDA dispatch.

Two entry points have a config but no demo script of their own. `CBIG_DELSSOME_finetuner.py` uses
`config/example_finetune_pFIC.ini`, which expects the predictor-training example to have run first
so it has a checkpoint to load. `CBIG_DELSSOME_CMAES_main_pMFM.py` uses `config/example_CMAES_pMFM.ini`
and runs on the same `input/` as the other CMA-ES demos — see
[../../group_level/README.md](../../group_level/README.md) §2.1 for the command.

## Layout

```
examples/group_level/
├── README.md                        <- this file
├── scripts/                         <- the four demo shell scripts
├── config/                          <- example .ini configs
├── input/                           <- demo SC / FC / FCD-CDF / myelin / gradient
├── reference_output/                <- reference snapshots for the checkers
├── CBIG_DELSSOME_check_*.py         <- one checker per demo
└── output/                          <- created at run time (gitignored)
```

To re-run cleanly, `rm -rf output/`. `input/` and `reference_output/` are read-only release
artefacts.
