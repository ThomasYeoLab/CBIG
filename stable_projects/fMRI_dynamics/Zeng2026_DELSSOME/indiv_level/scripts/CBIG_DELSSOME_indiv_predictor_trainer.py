# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
Step 2 — train the individual-level DELSSOME cost predictor.

Reads the per-subject ``param_save_epoch*.pth`` dumps produced by Step 1
(``CBIG_DELSSOME_indiv_cmaes.py --engine euler``) and trains a
:class:`MulModLossPredictorMfm` to predict the FC corr/L1/KS costs from
``(parameters, SC, FC, FCD-PDF)``.

Hyper-parameters are tuned with Optuna (``n_trials`` trials, capped at
``walltime``). The best checkpoint per trial is stored under
``LOG_DIR/all/train_indiv_DELSSOME-pfic/train/trial<N>/seed<seed>/lightning_logs/version_<k>/checkpoints/``.

A cost predictor is specific to a parcellation and FCD binning, so this is only
worth running for a large cohort on a new atlas -- see ``../README.md`` section 6
for the full three-step recipe. For fewer than ~20 subjects, skip the predictor
and run ``CBIG_DELSSOME_indiv_estimate_EI.py`` with pure Euler integration.

Usage (a small, fast configuration)::

    python -m scripts.CBIG_DELSSOME_indiv_predictor_trainer \\
        --trial 1 --seed 1 --ds-names MYSTUDY \\
        --model-type pfic-indiv-pred-multimodal \\
        --max-epochs 10 --n-trials 2 \\
        --walltime 00:30:00 --batch-size 256 --num-workers 0 \\
        --n-sample-per-phase train2_val2_test2

The paper's predictor used ``--max-epochs 15 --n-trials 50 --batch-size 2048
--n-sample-per-phase train40_val20_test20``, taking ~65 h on one RTX 3090.
"""

import argparse
from lightning import seed_everything

from DELSSOME_indiv.CBIG_file_utils import get_run_dir
from DELSSOME_indiv.CBIG_dl_train_utils import (
    walltime_to_seconds,
    config_gpu,
    tune,
    get_n_sample_per_phase,
    model_type_to_target,
)


def build_parser():
    p = argparse.ArgumentParser(description="Train the DELSSOME cost predictor.")
    p.add_argument("--trial", type=int, required=True)
    p.add_argument("--seed", type=int, required=True, help="Lightning + Optuna seed; also indexes the log folder.")
    p.add_argument("--model-type",
                   default="pfic-indiv-pred-multimodal",
                   help="'pfic-indiv-pred-multimodal' (default, transformer) "
                   "or 'pfic-indiv-pred' (MLP-only baseline).")
    p.add_argument("--max-epochs", type=int, default=15)
    p.add_argument("--n-trials", type=int, default=2, help="Number of Optuna trials.")
    p.add_argument("--walltime", default="120:00:00", help="Optuna study timeout, format HH:MM:SS.")
    p.add_argument("--batch-size", type=int, default=2048)
    p.add_argument("--num-workers", type=int, default=0)
    p.add_argument("--n-sample-per-phase",
                   default="train40_val20_test20",
                   help="Number of subjects per phase; supports stage chains "
                   "like 'train40_val20_test20-train40_val20'.")
    p.add_argument("--ds-names",
                   nargs="+",
                   default=None,
                   help="Datasets to train on. Defaults to ALL_TRAIN_DS "
                   "(12 datasets); pass e.g. '--ds-names HCP-YA' to restrict "
                   "to a single dataset when only that one was prepared.")
    p.add_argument("--enable-progress-bar", action="store_true")
    p.add_argument("--vis-study",
                   action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Write Optuna diagnostic plots (PNG export needs Chrome/kaleido). "
                   "Use --no-vis-study to skip (e.g. on machines without Chrome).")
    return p


def main():
    args = build_parser().parse_args()
    phase = "train"
    target = model_type_to_target(args.model_type)
    save_dir = get_run_dir("all", target, phase, args.trial, args.seed, None)

    seed_everything(args.seed)
    config_gpu()

    n_sample_per_phase = get_n_sample_per_phase(args.n_sample_per_phase)

    tune(
        model_type=args.model_type,
        max_epochs=args.max_epochs,
        n_trials=args.n_trials,
        timeout=walltime_to_seconds(args.walltime),
        save_dir=save_dir,
        seed=args.seed,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        n_sample_per_phase=n_sample_per_phase,
        ds_names=args.ds_names,
        enable_progress_bar=args.enable_progress_bar,
        vis_study=args.vis_study,
    )


if __name__ == "__main__":
    main()
