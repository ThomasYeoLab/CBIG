# DELSSOME: DEep Learning for Surrogate Statistics Optimization in MEan field modeling

## References

Zeng, T., Tian, F., Zhang, S., Li, X., Tan, A. P., Larsen, B., ... & Yeo, B. T. T. (2026). [Optimizing Biophysical Large-Scale Brain Circuit Models With Deep Neural Networks](https://www.biorxiv.org/content/10.1101/2025.04.07.647497v3). _Nature Methods_.

---

> **Stand-alone repository.** A self-contained copy of this project — including the `data/`
> training matrices and every file under `pretrained_models/`, which exceed this repository's
> 1 MB per-file limit and are therefore not committed here — is available at
> **https://github.com/TianCZeng/DELSSOME**. Clone that if you want everything in one place;
> otherwise see [Data](#data) below for what to download and where to put it.

Estimating the excitation/inhibition (E/I) ratio of the human cortex requires inverting a
biophysical circuit model — searching for the model parameters whose simulated fMRI best matches
real fMRI. The search is expensive because scoring a single candidate parameter set means
numerically integrating the model. DELSSOME replaces that integration with a transformer that
**predicts** the FC+FCD cost directly from the parameters, making each evaluation roughly 50×
cheaper.

This repository releases two pipelines. **The difference between them is the first thing to
understand:**

|                             | Group-level                                      | Individual-level                             |
| --------------------------- | ------------------------------------------------ | -------------------------------------------- |
| **What it estimates** | One E/I ratio map for a*group* of participants | One E/I ratio map**per participant**   |
| **What it fits to**   | **Group-averaged** FC and FCD CDF          | **Individual-specific** FC and FCD CDF |
| **Where**             | [group_level/](group_level/)                      | [indiv_level/](indiv_level/)                  |
| **Circuit models**    | pFIC, pMFM, Hopf                                 | pFIC                                         |
| **Extras**            | —                                               | GAMLSS normative lifespan trajectories       |

## Which part do I want?

| If you want to…                                                           | Go to                                                                                                                   |
| -------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| **Get an E/I ratio for your own subjects** (the most common request) | [examples/indiv_level/](examples/indiv_level/), then [indiv_level/README.md](indiv_level/README.md) §4                   |
| Fit a normative lifespan E/I trajectory to a cohort                        | [indiv_level/GAMLSS/](indiv_level/GAMLSS/)                                                                               |
| Invert a group-level model (pFIC, pMFM or Hopf), with or without DELSSOME  | [examples/group_level/](examples/group_level/), then [group_level/README.md](group_level/README.md)                       |
| Train your own DELSSOME cost predictor — rarely needed, since the release ships pretrained ones | [group_level/README.md](group_level/README.md) §4 (group) or [indiv_level/README.md](indiv_level/README.md) §7 (individual) |
| Reproduce the paper's full-dataset results (deferred to a future release)  | [replication/README.md](replication/README.md)                                                                           |

The single most common task has a single command. Given a functional-connectivity matrix and an
FCD cumulative distribution function for one subject:

```bash
cd indiv_level && export PYTHONPATH=$PWD
python -m scripts.CBIG_DELSSOME_indiv_estimate_EI \
    --fc my_data/sub-01_FC.csv --fcd-cdf my_data/sub-01_FCD_CDF.csv \
    --tr 2.0 --scan-length 600 \
    --out-dir results/sub-01
```

which writes `results/sub-01/EI_ratio.csv` — one row per cortical region. See
[indiv_level/README.md](indiv_level/README.md) for the input/output specification, and
[examples/indiv_level/](examples/indiv_level/) for a runnable version on bundled data.

## Overview

Instead of numerically integrating a candidate parameter set to find out how realistic the
resulting dynamics are, DELSSOME trains a deep neural network to predict the **FC+FCD cost** — the
surrogate statistics that quantify agreement between simulated and empirical fMRI — directly from
the parameters.

DELSSOME does not change the optimizer. The CMA-ES search is identical either way; what changes is
how each candidate parameter set is **evaluated** — a single forward pass of the cost predictor
instead of Euler integration of the circuit model. Following the paper, only the training stage
uses the predictor: validation and test always evaluate by Euler simulation, so the parameter set
you end up with has been verified by honest simulation. That also means the total wall-clock time
per subject is minutes to hours, not seconds.

### Supported circuit models

| Model          | Description                            | Group-level script                   | Individual-level |
| -------------- | -------------------------------------- | ------------------------------------ | ---------------- |
| **pFIC** | Parametric Feedback Inhibition Control | `CBIG_DELSSOME_CMAES_main_pFIC.py` | yes              |
| **pMFM** | Parametric Mean-Field Model            | `CBIG_DELSSOME_CMAES_main_pMFM.py` | no               |
| **Hopf** | Hopf Bifurcation Model                 | `CBIG_DELSSOME_CMAES_main_Hopf.py` | no               |

Optimization is CMA-ES (Covariance Matrix Adaptation Evolution Strategy), scoring candidates in
either `Euler` mode (direct simulation) or `DELSSOME` mode (the cost predictor).

All released models and pretrained predictors use **DK68** — the 68 cortical regions of the
Desikan–Killiany atlas (Desikan et al., 2006). A 100-region alternative is reported in the paper's
supplement but is not part of this release.

## Directory structure

Relative to the folder holding this README:

```
├── README.md                  # this file
├── group_level/               # group-level pipeline (pFIC, pMFM, Hopf)
│   ├── models/                # circuit models + DELSSOME predictor architecture
│   ├── process/               # data preprocessing and dataset classes
│   └── utils/                 # FC/FCD and parameter helpers
├── indiv_level/               # individual-level pipeline (pFIC)
│   ├── DELSSOME_indiv/        # Python package (model, predictor, CMA-ES, datasets)
│   ├── scripts/               # CLI entry points -- start with estimate_EI
│   ├── configs/               # INI configs
│   └── GAMLSS/                # normative lifespan trajectories (R)
├── examples/                  # runnable demos for BOTH pipelines  <- start here
│   ├── group_level/
│   └── indiv_level/
├── unit_tests/                # regression tests
├── replication/               # project configuration; full-dataset workflow deferred
├── data/                      # group-level predictor training data      [see Data]
└── pretrained_models/         # pretrained cost predictors               [see Data]
    ├── group_level/           # pFIC / pMFM / Hopf predictors (.pth)
    └── indiv_level/           # individual-level pFIC predictor (.ckpt)
```

## Data

This repository ships the code plus lightweight examples only. Files larger than this
repository's 1 MB per-file limit are **not committed here**:

| Not committed here | Needed by |
|---|---|
| `data/bootstrap_mats/`, `data/param_dataset/` | `CBIG_DELSSOME_predictor_train_example.sh` |
| `pretrained_models/group_level/` | `CBIG_DELSSOME_CMAES_pFIC_DL_example.sh` |
| `pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt` | the DELSSOME-accelerated part of `CBIG_DELSSOME_indiv_EI_example.sh` (optional) |

Take them from the stand-alone repository linked at the top of this page
(**https://github.com/TianCZeng/DELSSOME**) and place them at the same relative paths.

Everything else runs out of the box. In particular the individual-level example produces E/I
ratios with pure Euler integration, which needs no predictor at all — the download only adds the
accelerated variant. Inside CBIG, `unit_tests/` instead source their heavy data from
`$CBIG_TESTDATA_DIR/stable_projects/fMRI_dynamics/Zeng2026_DELSSOME/`
(see `unit_tests/README.md`).

### Example data provenance

**No empirical or otherwise access-controlled data is distributed in this repository.** Every
bundled example input is synthetic:

- `examples/indiv_level/input/` — the four `sub-000` … `sub-003` FC matrices, FCD CDFs, group SC,
  myelin and RSFC-gradient maps, and demographics are a deterministic function of
  `numpy.random.default_rng(42)`. Subject IDs are loop counters, not de-identified real
  identifiers.
- `examples/indiv_level/input/gamlss_input.csv` and `gamlss_new_cohort.csv` — synthesised by
  `indiv_level/GAMLSS/CBIG_generate_fake_gamlss_data.py`, with `fake-NNNN` IDs.
- `examples/group_level/input/` — group-level demo matrices, not traceable to any participant.

Because the inputs are synthetic, the numbers the examples produce are **not scientific results**.
They exist to show that the pipeline runs and to detect regressions.

Empirical data are governed by the original dataset providers' Data Use Agreements and are not
redistributed here. The pretrained predictors contain network weights and architecture
hyper-parameters only; `archive/CBIG_DELSSOME_sanitize_ckpt.py` is used before release to strip
optimizer state and the absolute filesystem paths that PyTorch Lightning embeds in checkpoints.

## Environment setup

For the **group-level** pipeline:

```bash
conda env create -f replication/config/CBIG_DELSSOME_python_env.yml
conda activate CBIG_DELSSOME
```

This installs Python 3.12, PyTorch 2.5.1 (CUDA 12.4), and all dependencies.

For the **individual-level** pipeline, use [indiv_level/environment.yml](indiv_level/environment.yml),
which additionally pulls in R + `r-gamlss` for the lifespan fit:

```bash
conda env create -f indiv_level/environment.yml
conda activate CBIG_DELSSOME_indiv
```

Either installation may take a few minutes.

## Quick start

```bash
# individual-level: your own FC / FCD-CDF -> per-region E/I ratio
bash examples/indiv_level/CBIG_DELSSOME_indiv_EI_example.sh

# group-level: CMA-ES inversion of the pFIC model
bash examples/group_level/scripts/CBIG_DELSSOME_CMAES_pFIC_example.sh
```

See [examples/README.md](examples/README.md) for all demos and their run times.

## Code Release

### Download stand-alone repository

Since the whole Github repository is too big, we provide a stand-alone version of only this project
and its dependencies. The stand-alone repository additionally carries the `data/` training matrices
and the `pretrained_models/` checkpoints, which exceed this repository's per-file size limit and are
therefore not committed here. To download this stand-alone repository, visit this link:
[https://github.com/TianCZeng/DELSSOME](https://github.com/TianCZeng/DELSSOME)

### Download whole repository

If you want to use the code from our lab's other stable projects (other than Zeng2026_DELSSOME),
you would want to download the whole CBIG repository.

- To download the version of the code that was last tested, you can either

    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.38.0-Zeng2026_DELSSOME](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.38.0-Zeng2026_DELSSOME)

    or

    - run the following command, if you have Git installed

    ```
    git checkout -b Zeng2026_DELSSOME v0.38.0-Zeng2026_DELSSOME
    ```

## Updates

- Release v0.38.0 (04/08/2026): Initial release of Zeng2026_DELSSOME

## License

This project is released under the MIT License, see
[https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md](https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md).

## Bugs and questions

Please contact Tianchu Zeng (tianchuzeng@gmail.com), Fang Tian (tianfang2001@gmail.com), or Prof. Thomas Yeo (yeoyeo02@gmail.com)
