# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
DELSSOME cost predictors + biophysical dynamic models.

* :class:`DELSSOME_indiv.models.CBIG_dynamic_model.pFIC` — the FIC mean-field
  model (Deco et al., 2014) with Euler integration and hemodynamic forward model.
* :class:`DELSSOME_indiv.models.CBIG_DELSSOME_predictor.MulModLossPredictorMfm`
  — the transformer cost predictor that replaces Euler integration inside CMA-ES.
"""

from DELSSOME_indiv.models.CBIG_dynamic_model import (
    DynamicModel,
    pFIC,
    AVAIL_MODELS,
    reshape_dup_sim_res,
    get_valid_params_after_dup,
)
from DELSSOME_indiv.models.CBIG_DELSSOME_predictor import (
    MLP,
    ParamEmbedder,
    SCEmbedder,
    FCEmbedder,
    FCDPDFEmbedder,
    LossPredictor,
    LossPredictorMfm,
    MulModLossPredictorMfm,
)

__all__ = [
    "DynamicModel",
    "pFIC",
    "AVAIL_MODELS",
    "reshape_dup_sim_res",
    "get_valid_params_after_dup",
    "MLP",
    "ParamEmbedder",
    "SCEmbedder",
    "FCEmbedder",
    "FCDPDFEmbedder",
    "LossPredictor",
    "LossPredictorMfm",
    "MulModLossPredictorMfm",
]
