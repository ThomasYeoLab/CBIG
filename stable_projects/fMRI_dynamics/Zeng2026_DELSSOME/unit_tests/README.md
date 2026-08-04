# Unit tests (internal use only)

Regression test for both Zeng2026_DELSSOME pipelines. It runs the bundled example scripts under
`../examples/` and asserts each completes with exit code 0. The numeric comparison against each
example's `reference_output/` is done by the examples' own checker scripts, which the example shell
scripts invoke themselves.

## Usage

```
sh CBIG_DELSSOME_unit_test.sh
```

or directly:

```
source activate CBIG_DELSSOME
python -u test_CBIG_DELSSOME_unit_test.py
```

## Tests

| Test | Example | Heavy data needed |
| --- | --- | --- |
| `test_cmaes_pfic` | `examples/group_level/scripts/CBIG_DELSSOME_CMAES_pFIC_example.sh` | no |
| `test_cmaes_hopf` | `examples/group_level/scripts/CBIG_DELSSOME_CMAES_Hopf_example.sh` | no |
| `test_cmaes_pfic_dl` | `examples/group_level/scripts/CBIG_DELSSOME_CMAES_pFIC_DL_example.sh` | `pretrained_models/` |
| `test_predictor_train` | `examples/group_level/scripts/CBIG_DELSSOME_predictor_train_example.sh` | `data/` |
| `test_indiv_estimate_ei` | `examples/indiv_level/CBIG_DELSSOME_indiv_EI_example.sh` | no |
| `test_indiv_gamlss` | `examples/indiv_level/CBIG_DELSSOME_indiv_gamlss_example.sh` | no |

## Prerequisites

- The `CBIG_DELSSOME` conda environment (`CBIG_DELSSOME_python_env.yml`, next to this file;
  identical to `../replication/config/CBIG_DELSSOME_python_env.yml`) for the four group-level
  tests.
- The `CBIG_DELSSOME_indiv` environment (`../indiv_level/environment.yml`) for the two
  individual-level tests. `test_indiv_gamlss` additionally wants `Rscript` with `gamlss`; without
  it the example exits 0 after printing a `[skip]` notice, so the test passes vacuously.
- Heavy data for two of the six tests, resolved in this order: (1) `../data/` and
  `../pretrained_models/` if present in the repo; (2) otherwise staged via symlink from
  `$CBIG_TESTDATA_DIR`; (3) otherwise the test is skipped with a message. Only symlinks created by
  the test are removed afterwards, so real in-repo directories are never clobbered.

`test_indiv_estimate_ei` runs pure Euler and therefore needs no predictor. The example adds a
DELSSOME-accelerated run only when `pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt` exists, and
skips that part cleanly otherwise — so the test covers the same code path either way.

`../data/` and `../pretrained_models/` are **not** committed in this repository. Place the data
under `$CBIG_TESTDATA_DIR` mirroring the in-repo layout:

```
$CBIG_TESTDATA_DIR/stable_projects/fMRI_dynamics/Zeng2026_DELSSOME/
├── data/
│   ├── bootstrap_mats/{train,val,test}.mat
│   └── param_dataset/{train,val,test}/{group0,group1}/...
└── pretrained_models/
    ├── group_level/pFIC_DK68_predictor.pth
    └── indiv_level/pFIC_DK68_predictor.ckpt  (optional; extends test_indiv_estimate_ei)
```

Alternatively, obtain them from the stand-alone public repository
(**https://github.com/TianCZeng/DELSSOME**) and place them at `../data/` and
`../pretrained_models/`.

Runtime output is written under each example's own `output/` directory (git-ignored). These tests
cover the bundled examples only; reproducing the paper's full-dataset results is deferred to a
future release — see [../replication/README.md](../replication/README.md).
