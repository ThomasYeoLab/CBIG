# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Project-wide constants and default paths.

The project directory is taken from the ``DELSSOME_PROJECT_DIR`` environment
variable, falling back to the parent directory of this package (i.e. the
``indiv_level/`` folder) so the bundled examples work with no configuration.

All downstream code computes its inputs/outputs relative to ``PROJECT_DIR``,
``DATA_DIR``, ``LOG_DIR`` and ``CONFIG_DIR``, so re-routing those four
variables is enough to make the pipeline run anywhere.
"""

import os
import torch

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
# ``indiv_level/`` directory (parent of ``DELSSOME_indiv/``).
_DEFAULT_PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Top-level project directory used for resolving data/logs/configs.
#: Override by exporting ``DELSSOME_PROJECT_DIR`` before importing the package.
PROJECT_DIR = os.environ.get("DELSSOME_PROJECT_DIR", _DEFAULT_PROJECT_DIR)

CONFIG_DIR = os.environ.get("DELSSOME_CONFIG_DIR", os.path.join(PROJECT_DIR, "configs"))
DATA_DIR = os.environ.get("DELSSOME_DATA_DIR", os.path.join(PROJECT_DIR, "examples", "data"))
LOG_DIR = os.environ.get("DELSSOME_LOG_DIR", os.path.join(PROJECT_DIR, "examples", "logs"))
TMP_DIR = os.environ.get("DELSSOME_TMP_DIR", os.path.join(PROJECT_DIR, "examples", "tmp"))

# Where group-level HCP-YA reference statistics (SC, myelin, RSFC gradient) live.
# In the example data shipped with this release, these come from
# ``examples/data/HCP-YA/pFIC_input/<atlas>/``.
HCPYA_pFIC_INPUT = os.path.join(DATA_DIR, "HCP-YA", "pFIC_input")

# --------------------------------------------------------------------------- #
# Parcellation
# --------------------------------------------------------------------------- #
NUM_ROI = 68  # Desikan-Killiany cortical parcellation, minus 4 ignored ROIs
ATLAS_TO_N_ROI = {"DK68": 68}

# --------------------------------------------------------------------------- #
# Datasets
# --------------------------------------------------------------------------- #
# Datasets used to TRAIN the individual-level DELSSOME cost predictor in the
# paper (Section 2.6). Listed here so that ``IndivMulModPredDs2014`` can
# attach a numeric dataset index to every (parameter, loss) training sample.
ALL_TRAIN_DS = [
    "ABCD", "GUSTO", "devCCNP", "PNC", "HCP-D", "HCP-YA", "HCP-A", "eNKI",
    "SINGER", "SG70", "MACC", "ADNI",
]
ALL_EVAL_DS = ["LIFE", "TCP", "UKB"]
ALL_DS = ALL_TRAIN_DS + ALL_EVAL_DS
DS_IDX = {ds: idx for idx, ds in enumerate(ALL_DS)}

#: Index handed to datasets that are not one of the paper's 15. The predictor
#: only uses this index for logging (it is not an input to ``forward``), so
#: training on a cohort of your own is unaffected by not being listed above.
OTHER_DS_IDX = len(ALL_DS)


def get_ds_idx(ds_name):
    """Numeric dataset index, falling back to :data:`OTHER_DS_IDX` with a warning."""
    if ds_name in DS_IDX:
        return DS_IDX[ds_name]
    print(f"[warn] '{ds_name}' is not one of the {len(ALL_DS)} datasets used in the "
          f"paper; assigning the generic index {OTHER_DS_IDX}. This only affects "
          "per-dataset logging, not the trained predictor.")
    return OTHER_DS_IDX


# --------------------------------------------------------------------------- #
# Optimisation phases
# --------------------------------------------------------------------------- #
#: Maps a current phase to the previous phase whose best parameters we
#: bootstrap from. ``train -> train`` because the first phase has no parent.
PREV_PHASE = {"train": "train", "val": "train", "test": "val"}

# --------------------------------------------------------------------------- #
# Numerical
# --------------------------------------------------------------------------- #
#: Default torch dtype used throughout. Float64 to keep the Euler integration
#: numerically stable; cast to float32 for the deep-learning training pipeline.
DEFAULT_DTYPE = torch.float64

#: Loss names (in canonical order) used to construct the FC+FCD cost.
ALL_LOSSES = ["total_loss", "corr_loss", "l1_loss", "ks_loss"]
