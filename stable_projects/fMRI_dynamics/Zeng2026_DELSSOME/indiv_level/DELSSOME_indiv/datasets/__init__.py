# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
DataModules for training the DELSSOME cost predictor.

The training data are produced by per-subject Euler CMA-ES runs
(``CBIG_DELSSOME_indiv_cmaes.py --engine euler``), which drop one
``param_save_epoch<k>.pth``
file per training epoch under ``LOG_DIR/<ds>/indiv_Euler-pfic/train/...``.
These files contain the (parameter, valid mask, three losses) tuples that the
DELSSOME predictor learns to regress.
"""

from DELSSOME_indiv.datasets.CBIG_dl_dataset import (
    n_sample_per_phase_to_list,
    n_sample_per_phase_to_file_name,
    ParamDataset,
    DataModule,
    PredictorDataset,
    MulModPredDs,
    IndivMulModPredDs2014,
    IndivMulModPredDm2014,
    IndivPredDataset,
    IndivPredDataModule,
)

__all__ = [
    "n_sample_per_phase_to_list",
    "n_sample_per_phase_to_file_name",
    "ParamDataset",
    "DataModule",
    "PredictorDataset",
    "MulModPredDs",
    "IndivMulModPredDs2014",
    "IndivMulModPredDm2014",
    "IndivPredDataset",
    "IndivPredDataModule",
]
