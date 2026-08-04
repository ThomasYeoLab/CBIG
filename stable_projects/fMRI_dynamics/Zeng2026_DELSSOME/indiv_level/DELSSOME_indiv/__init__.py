# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
DELSSOME — DEep Learning for Surrogate Statistics Optimization in MEan field modeling.

This package contains the *individual-level* DELSSOME implementation released
alongside Zeng et al. (2026). It provides:

* A PyTorch implementation of the feedback-inhibition-control (FIC) mean-field
  model of Deco et al. (2014) (`DELSSOME_indiv.models.CBIG_dynamic_model.pFIC`).
* The transformer-based DELSSOME cost predictor that replaces costly Euler
  integration during CMA-ES parameter inversion
  (`DELSSOME_indiv.models.CBIG_DELSSOME_predictor.MulModLossPredictorMfm`).
* The CMA-ES driver that calls either Euler integration or DELSSOME to evaluate
  candidate parameter sets (`DELSSOME_indiv.optimizers.CBIG_cmaes`).
* A PyTorch-Lightning data pipeline for assembling DELSSOME training data from
  per-subject Euler CMA-ES runs (`DELSSOME_indiv.datasets.CBIG_dl_dataset`).
* Utility helpers (paths, FC/FCD, losses, etc.).

The package is intentionally a near-verbatim refactor of the research codebase
to keep behaviour identical to what is reported in the paper, while removing
hard-coded cluster paths so that it can be exercised on a local machine.
"""

from DELSSOME_indiv.CBIG_constants import (
    PROJECT_DIR,
    CONFIG_DIR,
    DATA_DIR,
    LOG_DIR,
    NUM_ROI,
    ATLAS_TO_N_ROI,
    DEFAULT_DTYPE,
    PREV_PHASE,
)

__all__ = [
    "PROJECT_DIR",
    "CONFIG_DIR",
    "DATA_DIR",
    "LOG_DIR",
    "NUM_ROI",
    "ATLAS_TO_N_ROI",
    "DEFAULT_DTYPE",
    "PREV_PHASE",
]
