# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Step 2b (optional) — evaluate a trained DELSSOME checkpoint on the test split.

This loads a Lightning checkpoint produced by
``CBIG_DELSSOME_indiv_predictor_trainer.py`` and runs ``Trainer.test`` against
the ``IndivMulModPredDs2014`` test dataset. The per-prediction CSV
(``test_pred.csv``) is dropped next to the checkpoint.

Example::

    python -m scripts.CBIG_DELSSOME_indiv_predictor_tester \\
        --ckpt-path <run_dir>/lightning_logs/version_0/checkpoints/last.ckpt \\
        --model-type pfic-indiv-pred-multimodal \\
        --n-sample-per-phase train2_val2_test2
"""

import argparse
import os
from DELSSOME_indiv.CBIG_dl_train_utils import (
    test,
    get_n_sample_per_phase,
    config_gpu,
    model_type_to_target,
)
from DELSSOME_indiv.CBIG_file_utils import get_run_dir


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--ckpt-path",
        required=True,
        help="Path to a Lightning .ckpt produced by "
        "CBIG_DELSSOME_indiv_predictor_trainer.py.")
    p.add_argument("--model-type", default="pfic-indiv-pred-multimodal")
    p.add_argument("--trial", type=int, default=1)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--batch-size", type=int, default=2048)
    p.add_argument("--num-workers", type=int, default=0)
    p.add_argument("--n-sample-per-phase", default="train40_val20_test20")
    p.add_argument("--ds-names",
                   nargs="+",
                   default=None,
                   help="Datasets to evaluate on. Defaults to ALL_TRAIN_DS.")
    p.add_argument("--enable-progress-bar", action="store_true")
    return p


def main():
    args = build_parser().parse_args()
    target = model_type_to_target(args.model_type)
    save_dir = get_run_dir("all", target, "test", args.trial, args.seed, None)
    os.makedirs(save_dir, exist_ok=True)
    config_gpu()
    test(
        model=args.ckpt_path,
        model_type=args.model_type,
        save_dir=save_dir,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        n_sample_per_phase=get_n_sample_per_phase(args.n_sample_per_phase),
        ds_names=args.ds_names,
        enable_progress_bar=args.enable_progress_bar,
    )


if __name__ == "__main__":
    main()
