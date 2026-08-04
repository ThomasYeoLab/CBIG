# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Per-subject CMA-ES for the individual-level pFIC model.

One script, two engines:

* ``--engine euler`` (default) scores every candidate parameter set by Euler
  integration. This is the classical approach; it needs no trained predictor.
  Running it over a pilot set of subjects is also how you GENERATE the training
  data for the DELSSOME cost predictor.
* ``--engine delssome`` replaces the Euler inner loop of the *training* stage
  with a forward pass through a trained DELSSOME predictor, giving the speed-up
  reported in the paper. Validation and test keep Euler integration, as the
  paper recommends, so they are unaffected.

Three stages run per subject, chained through files on disk:

* ``train`` — ``num_epochs`` epochs of CMA-ES, writing ``param_save_epoch{k}.pth``
  and then ``best_from_train.pth``.
* ``val``   — re-simulate ``best_from_train.pth`` with the validation dt
  (writes ``best_from_val.pth``).
* ``test``  — re-simulate the best validation parameter with the test dt
  (writes ``best_from_test.pth``, which holds the E/I-ratio estimates).

Use ``--mode all`` to run the three in sequence.

.. note::
   ``--mode train`` starts from a clean slate: it deletes the run directory
   before epoch 0. Pass ``--next-epoch`` to resume instead.

Usage::

    # Euler, all three stages, one subject
    python -m scripts.CBIG_DELSSOME_indiv_cmaes \\
        --ds-name HCP-YA --sub-id sub-000 --mode all --num-epochs 100

    # DELSSOME-accelerated training, Euler validation/test
    python -m scripts.CBIG_DELSSOME_indiv_cmaes \\
        --ds-name HCP-YA --sub-id sub-000 --mode all --num-epochs 100 \\
        --engine delssome \\
        --delssome-ckpt ../pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt

If you just want an E/I ratio from your own FC / FCD-CDF CSVs and do not care
about the log-directory layout, use ``CBIG_DELSSOME_indiv_estimate_EI.py``
instead — it wraps this module behind a much smaller set of flags.
"""

import argparse
import datetime
import sys
import numpy as np
import torch

from DELSSOME_indiv.CBIG_misc_utils import set_torch_default
from DELSSOME_indiv.CBIG_file_utils import (
    get_stats,
    get_config,
    get_prev_phase_best_params_path,
    get_curr_phase_save_dir,
    get_best_train_params,
)
from DELSSOME_indiv.optimizers import (
    ModelValidator,
    ModelTester,
    train_helper,
)

# --------------------------------------------------------------------------- #
# Phase routing
# --------------------------------------------------------------------------- #
VALID_MODES = (
    "all",  # train -> val -> test in sequence
    "train",  # CMA-ES, writes per-epoch dumps + best_from_train.pth
    "val",  # re-simulate best_from_train.pth, writes best_from_val.pth
    "test",  # re-simulate best_from_val.pth, writes best_from_test.pth
    "best_from_train",  # only rebuild best_from_train.pth from existing per-epoch dumps
)

#: Log-folder target per engine, so Euler and DELSSOME runs never overwrite
#: each other even with the same trial/seed indices.
ENGINE_TO_TARGET = {
    "euler": "indiv_Euler-pfic",
    "delssome": "apply_indiv_DELSSOME-pfic",
}


def extract_phase(mode):
    if mode in ("train", "val", "test"):
        return mode
    if mode == "best_from_train":
        return "train"
    raise ValueError(f"Invalid mode '{mode}'.")


def get_dynamic_model_name_from_target(target):
    return "mfm2013" if "2013" in target else "pfic"


# --------------------------------------------------------------------------- #
# Main driver
# --------------------------------------------------------------------------- #
def apply(ds_name,
          target,
          mode,
          trial_idx,
          seed_idx,
          group_or_sub_id,
          group_type=None,
          atlas="DK68",
          num_epochs=100,
          seed=None,
          agg_seeds_num=None,
          use_DELSSOME=False,
          delssome_ckpt=None,
          next_epoch=0,
          config=None,
          emp_stats=None):
    """Run one stage of CMA-ES for one subject.

    Args:
        ds_name: dataset name (matches the folders under ``DATA_DIR``).
        target:  log target — see :data:`ENGINE_TO_TARGET`.
        mode:    one of :data:`VALID_MODES`, excluding ``all``.
        trial_idx / seed_idx: arbitrary integers used to organise runs.
        group_or_sub_id:      subject identifier (string).
        use_DELSSOME:         if True, the Euler step of the training stage is
                              replaced with the DELSSOME cost predictor
                              (requires ``delssome_ckpt``).
        delssome_ckpt:        path to a trained DELSSOME Lightning checkpoint.
        next_epoch:           resume the training stage from this epoch instead
                              of wiping the run directory and starting at 0.
        config:               pre-built ``ConfigParser``. When None it is
                              resolved from ``ds_name``/``atlas`` on disk.
        emp_stats:            pre-built empirical-statistics dict. When None it
                              is loaded from ``DATA_DIR`` by convention. Both
                              hooks exist so ``CBIG_DELSSOME_indiv_estimate_EI``
                              can supply user CSVs without a data tree.

    Returns:
        int: 0 on success; 1 or 2 when CMA-ES could not be initialised or broke
        mid-run (see :func:`DELSSOME_indiv.optimizers.train_helper`).
    """
    print("Input arguments:")
    print(f"- ds_name: {ds_name}, target: {target}, mode: {mode}, "
          f"trial_idx: {trial_idx}, seed_idx: {seed_idx}, "
          f"group_or_sub_id: {group_or_sub_id}, group_type: {group_type}, "
          f"atlas: {atlas}, num_epochs: {num_epochs}, "
          f"seed: {seed}, agg_seeds_num: {agg_seeds_num}, "
          f"use_DELSSOME: {use_DELSSOME}")
    set_torch_default()

    dynamic_model_name = get_dynamic_model_name_from_target(target)
    print(f"Using dynamic model: {dynamic_model_name}")
    if config is None:
        config = get_config(ds_name,
                            atlas,
                            group_or_sub_id=group_or_sub_id,
                            group_type=group_type,
                            dynamic_model_name=dynamic_model_name)
    phase = extract_phase(mode)
    trial_idx = int(trial_idx)
    seed_idx = int(seed_idx)

    prev_phase_best_params_path = get_prev_phase_best_params_path(
        ds_name=ds_name,
        target=target,
        curr_phase=phase,
        trial_idx=trial_idx,
        seed_idx=seed_idx,
        curr_group_or_sub_id=group_or_sub_id,
        group_type=group_type,
        agg_seeds_num=agg_seeds_num,
    )
    curr_phase_save_dir = get_curr_phase_save_dir(
        ds_name=ds_name,
        target=target,
        phase=phase,
        trial_idx=trial_idx,
        seed_idx=seed_idx,
        group_or_sub_id=group_or_sub_id,
        agg_seeds_num=agg_seeds_num,
    )

    with_dl_stats = phase == "train"
    use_harmonized = group_type is not None
    if emp_stats is None:
        emp_stats = get_stats(ds_name=ds_name,
                              group_or_sub_id=group_or_sub_id,
                              group_type=group_type,
                              atlas=atlas,
                              with_myelin_gradient=True,
                              with_dl_stats=with_dl_stats,
                              use_harmonized=use_harmonized)

    # Tell the dynamic model where to find the DELSSOME checkpoint (if any).
    if use_DELSSOME and delssome_ckpt is not None:
        emp_stats["delssome_ckpt"] = delssome_ckpt

    if mode == "train":
        if use_DELSSOME:
            dl_epoch_range = np.arange(0, num_epochs)
            euler_epoch_range = []
        else:
            dl_epoch_range = []
            euler_epoch_range = np.arange(0, num_epochs)
        state = train_helper(
            config=config,
            emp_stats=emp_stats,
            train_save_dir=curr_phase_save_dir,
            num_epochs=num_epochs,
            dl_epoch_range=dl_epoch_range,
            euler_epoch_range=euler_epoch_range,
            opportunities=10,
            next_epoch=next_epoch,
            seed=seed,
        )
        print("Exit state: ", state)
        if state != 0:
            return state
        get_best_train_params(curr_phase_save_dir,
                              num_epochs=num_epochs,
                              top_k_for_each_epoch=1)
        return state

    if mode == "best_from_train":
        get_best_train_params(curr_phase_save_dir,
                              num_epochs=num_epochs,
                              top_k_for_each_epoch=1)
        return 0

    if mode == "val":
        model_handler = ModelValidator(config, emp_stats,
                                       prev_phase_best_params_path,
                                       curr_phase_save_dir)
        model_handler.validate(seed=seed)
        return 0

    if mode == "test":
        model_handler = ModelTester(config, emp_stats,
                                    prev_phase_best_params_path,
                                    curr_phase_save_dir)
        model_handler.test(seed=seed)
        return 0

    raise ValueError(f"Mode '{mode}' not implemented.")


def apply_all_stages(**kwargs):
    """Run ``train`` -> ``val`` -> ``test`` for one subject, stopping on failure.

    ``emp_stats`` is reloaded per stage because only the training stage needs the
    vectorised DL fields, and the validation/test stages use a different dt.
    """
    stage_emp_stats = kwargs.pop("emp_stats_per_stage", None)
    for stage in ("train", "val", "test"):
        print(f"\n=== stage: {stage} ===", flush=True)
        stage_kwargs = dict(kwargs, mode=stage)
        if stage_emp_stats is not None:
            stage_kwargs["emp_stats"] = stage_emp_stats[stage]
        state = apply(**stage_kwargs)
        if state != 0:
            print(f"[error] stage '{stage}' failed with state {state}; stopping.")
            return state
    return 0


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Per-subject CMA-ES for the individual-level pFIC model "
        "(Euler or DELSSOME-accelerated).")
    p.add_argument("--ds-name",
                   required=True,
                   help="Dataset name (folder name under DATA_DIR).")
    p.add_argument("--sub-id", required=True, help="Subject ID.")
    p.add_argument("--engine",
                   default="euler",
                   choices=sorted(ENGINE_TO_TARGET),
                   help="Cost evaluation inside the CMA-ES training stage. "
                   "'euler' needs no predictor; 'delssome' needs "
                   "--delssome-ckpt (default: euler).")
    p.add_argument("--mode",
                   default="all",
                   choices=VALID_MODES,
                   help="Stage to run, or 'all' for train -> val -> test "
                   "(default: all).")
    p.add_argument("--target",
                   default=None,
                   help="Log-folder target. Defaults to the value implied by "
                   "--engine, so Euler and DELSSOME runs stay separate.")
    p.add_argument("--trial", default=1, type=int,
                   help="Folder label only; groups runs under trial<N>.")
    p.add_argument("--seed-idx", default=1, type=int,
                   help="Folder label only; groups runs under seed<N>. "
                   "Not the RNG seed -- that is --seed.")
    p.add_argument("--group-type",
                   default=None,
                   help="None for individual-level (default).")
    p.add_argument("--atlas", default="DK68")
    p.add_argument("--num-epochs", default=100, type=int)
    p.add_argument("--next-epoch", default=0, type=int,
                   help="Resume the training stage from this epoch. The "
                   "default (0) wipes the run directory first.")
    p.add_argument("--seed",
                   default=None,
                   type=int,
                   help="RNG seed (None = sample one).")
    p.add_argument("--agg-seeds-num", default=None, type=int,
                   help="Merge this many consecutive training seeds into one "
                   "parameter pool at validation/test time.")
    p.add_argument(
        "--delssome-ckpt",
        default=None,
        help="Path to a trained DELSSOME .ckpt (required for --engine delssome).")
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    use_delssome = args.engine == "delssome"
    target = args.target or ENGINE_TO_TARGET[args.engine]
    if use_delssome and args.mode in ("all", "train") and not args.delssome_ckpt:
        parser.error(
            "--engine delssome needs --delssome-ckpt <path-to-predictor.ckpt>. "
            "Use --engine euler to run without a trained predictor.")

    print(datetime.datetime.now(), ": MAIN PROGRAM START", flush=True)
    if torch.cuda.is_available():
        print("Current GPU:",
              torch.cuda.get_device_name(torch.cuda.current_device()))

    runner = apply_all_stages if args.mode == "all" else apply
    kwargs = dict(
        ds_name=args.ds_name,
        target=target,
        trial_idx=args.trial,
        seed_idx=args.seed_idx,
        group_or_sub_id=args.sub_id,
        group_type=args.group_type,
        atlas=args.atlas,
        num_epochs=args.num_epochs,
        seed=args.seed,
        agg_seeds_num=args.agg_seeds_num,
        use_DELSSOME=use_delssome,
        delssome_ckpt=args.delssome_ckpt,
        next_epoch=args.next_epoch,
    )
    if args.mode != "all":
        kwargs["mode"] = args.mode
    # Non-zero exit on failure: the old scripts printed the state and exited 0,
    # so a broken CMA-ES run looked like a success to any calling shell script.
    sys.exit(runner(**kwargs))


if __name__ == "__main__":
    main()
