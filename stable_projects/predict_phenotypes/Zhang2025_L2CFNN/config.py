#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
"""
Written by Chen Zhang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
from pathlib import Path


class global_config:
    # Determine root path: Check environment variable first, then default.
    _default_root = Path(
        __file__
    ).parent / "unit_tests" # Directory containing this config.py
    root_path = Path(os.getenv("CBIG_L2C_ROOT_PATH", _default_root)).resolve()

    # Check if we are running inside the unit test environment
    is_unit_test = "unit_tests" in str(root_path)

    # Choose base path for data depending on test mode
    base_data_path = root_path / "output" if is_unit_test else root_path

    # Define paths relative to base_data_path
    raw_data_path = base_data_path / "raw_data"
    split_path = base_data_path / "data" / "folds_split"
    baseline_path = base_data_path / "data" / "baselines"
    fold_gen_path = base_data_path / "data" / "fold_gen"
    frog_path = base_data_path / "data" / "frog_feature"
    processed_frog = base_data_path / "data" / "processed_frog"
    clean_gt_path = base_data_path / "data" / "folds_split" / "clean_gt"
    minrnn_path = base_data_path / "data" / "minrnn_data"

    checkpoint_dir = base_data_path / "checkpoints"
    prediction_dir = base_data_path / "predictions"
    reference_dir = root_path / "ref_results"

    # Only for predict_only
    pretrained_dir = root_path / "pretrained_model"

    feature_dict = {
        "FULL": [
            "ADAS13",
            "MMSE",
            "CDR",
            "FDG",
            "ADAS11",
            "RAVLT_immediate",
            "RAVLT_learning",
            "RAVLT_perc_forgetting",
            "Ventricles",
            "Fusiform",
            "WholeBrain",
            "Hippocampus",
            "MidTemp",
            "ICV",
        ],
        "PARTIAL": [
            "MMSE",
            "CDR",
            "Ventricles",
            "Fusiform",
            "WholeBrain",
            "Hippocampus",
            "MidTemp",
            "ICV",
        ],
        "MinFull": [
            "CDRSB",
            "ADAS11",
            "ADAS13",
            "MMSE",
            "RAVLT_immediate",
            "RAVLT_learning",
            "RAVLT_forgetting",
            "RAVLT_perc_forgetting",
            "MOCA",
            "FAQ",
            "Entorhinal",
            "Fusiform",
            "Hippocampus",
            "ICV",
            "MidTemp",
            "Ventricles",
            "WholeBrain",
            "AV45",
            "FDG",
            "ABETA_UPENNBIOMK9_04_19_17",
            "TAU_UPENNBIOMK9_04_19_17",
            "PTAU_UPENNBIOMK9_04_19_17",
        ],
    }
