# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Estimate a per-region excitation/inhibition (E/I) ratio for ONE subject from
your own FC and FCD-CDF matrices.

This is the entry point to use if you have empirical data and want E/I numbers.
It hides the log-directory layout, the three-stage chaining, the subject-list
files and the dataset-name conventions that the replication scripts need, and
it writes plain CSVs instead of PyTorch pickles.

Example::

    cd indiv_level
    export PYTHONPATH=$PWD
    python -m scripts.CBIG_DELSSOME_indiv_estimate_EI \\
        --fc          my_data/sub-01_FC.csv \\
        --fcd-cdf     my_data/sub-01_FCD_CDF.csv \\
        --tr          2.0 \\
        --scan-length 600 \\
        --out-dir     results/sub-01

What happens internally: CMA-ES searches the pFIC parameter space for the
parameters whose simulated FC and FCD best match yours, over three stages --
``train`` (the search), ``val`` (re-score the per-epoch winners), and ``test``
(re-simulate the single best parameter set at a finer Euler step, and read off
the E/I ratio).

Pass ``--delssome-ckpt`` to run the *training* stage through the trained
DELSSOME cost predictor instead of Euler integration. Validation and test keep
Euler integration, following the paper. That means DELSSOME shortens the search
but not the final re-simulation, so expect minutes-to-hours per subject either
way -- it is not an interactive tool.

Inputs (all header-less CSV, N = number of ROIs):

    --fc          N x N     required   symmetric, diagonal 1
    --fcd-cdf     bins x 1  required   non-decreasing, last value 1
    --sc          N x N     optional   defaults to the bundled HCP-YA group SC
    --myelin      N x 1     optional   defaults to the bundled HCP-YA map
    --gradient    N x 1     optional   defaults to the bundled HCP-YA map

Outputs written to ``--out-dir``:

    EI_ratio.csv         region_index, EI_ratio, S_E, S_I, r_E
    best_parameters.csv  region_index, wEE, wEI, sigma, G
    fit_quality.csv      stage, corr_loss, l1_loss, ks_loss, total_loss
    run_config.ini       the fully resolved config, for provenance
    .work/               intermediate .pth files (safe to delete)
"""

import argparse
import datetime
import os
import shutil
import sys
import numpy as np
import pandas as pd
import torch

from DELSSOME_indiv.CBIG_constants import CONFIG_DIR
from DELSSOME_indiv.CBIG_misc_utils import set_torch_default
from DELSSOME_indiv.CBIG_file_utils import (
    build_config,
    build_stats_from_paths,
    get_best_params_file_path,
    get_best_train_params,
)
from DELSSOME_indiv.optimizers import ModelTester, ModelValidator, train_helper

#: Stage order. ``train`` searches, ``val`` re-scores, ``test`` re-simulates the
#: winner at the finer test dt and produces the E/I ratio.
STAGES = ("train", "val", "test")


# --------------------------------------------------------------------------- #
# Input validation
# --------------------------------------------------------------------------- #
def read_csv_matrix(path, flag):
    """Read a header-less numeric CSV, with an actionable error message."""
    if not os.path.exists(path):
        raise SystemExit(f"[error] {flag}: file not found: {path}")
    try:
        values = pd.read_csv(path, header=None).values.astype(float)
    except ValueError as exc:
        raise SystemExit(
            f"[error] {flag}: {path} is not a plain numeric CSV ({exc}). "
            "The file must have no header row and no index column.")
    return values


def validate_inputs(args):
    """Check every input up front and return the ROI count.

    Every failure here is something the user can fix from the message alone;
    the alternative is a shape error thousands of Euler steps later.
    """
    fc = read_csv_matrix(args.fc, "--fc")
    if fc.ndim != 2 or fc.shape[0] != fc.shape[1]:
        raise SystemExit(
            f"[error] --fc: expected a square N x N matrix, got {fc.shape}. "
            "Header rows or an index column are the usual cause.")
    n_roi = fc.shape[0]
    if not np.allclose(fc, fc.T, atol=1e-6):
        raise SystemExit("[error] --fc: matrix is not symmetric. Pass a "
                         "correlation matrix, not a time series.")

    fcd = read_csv_matrix(args.fcd_cdf, "--fcd-cdf")
    fcd = fcd.reshape(-1)
    if np.any(np.diff(fcd) < -1e-9):
        raise SystemExit("[error] --fcd-cdf: values must be non-decreasing. "
                         "This file should be a cumulative distribution "
                         "function, not a histogram or a PDF.")
    if not np.isclose(fcd[-1], 1.0, atol=1e-6):
        print(f"[warn] --fcd-cdf: last value is {fcd[-1]:.6g}, not 1.0; it "
              "will be rescaled so the CDF ends at 1.")

    for path, flag in ((args.sc, "--sc"), ):
        if path is None:
            continue
        mat = read_csv_matrix(path, flag)
        if mat.shape != (n_roi, n_roi):
            raise SystemExit(
                f"[error] {flag}: shape {mat.shape} does not match --fc "
                f"({n_roi} x {n_roi}).")

    for path, flag in ((args.myelin, "--myelin"), (args.gradient, "--gradient")):
        if path is None:
            continue
        vec = read_csv_matrix(path, flag).reshape(-1)
        if vec.size != n_roi:
            raise SystemExit(
                f"[error] {flag}: has {vec.size} values but --fc implies "
                f"{n_roi} ROIs.")

    if (args.myelin is None) != (args.gradient is None):
        raise SystemExit(
            "[error] --myelin and --gradient parameterise wEE/wEI/sigma "
            "together; pass both or neither.")
    if args.myelin is None and n_roi != 68:
        raise SystemExit(
            f"[error] --fc implies {n_roi} ROIs, but the bundled default myelin "
            "and RSFC-gradient maps are DK68 (68 regions). Pass --myelin and "
            "--gradient for your own parcellation.")

    return n_roi, fcd.size


def validate_ckpt(ckpt_path, n_roi, n_bins):
    """Check a DELSSOME checkpoint matches the data before we rely on it.

    The checkpoint records its own architecture, so this is cheap and turns a
    confusing failure deep inside the predictor's forward pass into one line.
    """
    if not os.path.exists(ckpt_path):
        raise SystemExit(f"[error] --delssome-ckpt: file not found: {ckpt_path}")
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    hparams = ckpt.get("hyper_parameters", {}).get("model_hparams", {})
    ckpt_n_roi = hparams.get("n_regions")
    ckpt_bins = hparams.get("bins")
    if ckpt_n_roi is not None and int(ckpt_n_roi) != n_roi:
        raise SystemExit(
            f"[error] --delssome-ckpt was trained for {ckpt_n_roi} ROIs but "
            f"your FC has {n_roi}. A predictor cannot be reused across "
            "parcellations -- train one on your own data, or drop "
            "--delssome-ckpt to run pure Euler.")
    if ckpt_bins is not None and int(ckpt_bins) != n_bins:
        raise SystemExit(
            f"[error] --delssome-ckpt expects {ckpt_bins} FCD-CDF bins but "
            f"your --fcd-cdf has {n_bins} rows.")
    print(f"DELSSOME checkpoint OK: {ckpt_n_roi} ROIs, {ckpt_bins} FCD bins")


# --------------------------------------------------------------------------- #
# Output writers
# --------------------------------------------------------------------------- #
def as_column(tensor):
    """Squeeze a ``[N, 1]`` or ``[N]`` tensor to a 1D numpy array."""
    return torch.as_tensor(tensor).detach().cpu().reshape(-1).numpy()


def write_ei_csv(best_from_test, out_dir):
    """Write the headline result: one row per region."""
    columns = {
        "EI_ratio": "EI_for_valid_params",
        "S_E": "S_E_for_valid_params",
        "S_I": "S_I_for_valid_params",
        "r_E": "r_E_for_valid_params",
    }
    data = {}
    for out_name, key in columns.items():
        if key in best_from_test:
            data[out_name] = as_column(best_from_test[key])
    if "EI_ratio" not in data:
        raise SystemExit(
            "[error] the test stage produced no E/I ratio. This usually means "
            "no candidate parameter set stayed in range -- try more epochs.")
    n_roi = len(data["EI_ratio"])
    frame = pd.DataFrame({"region_index": np.arange(n_roi), **data})
    path = os.path.join(out_dir, "EI_ratio.csv")
    frame.to_csv(path, index=False)
    print(f"Wrote {path}")
    return frame


def write_parameters_csv(best_from_test, n_roi, out_dir):
    """Split the 3N+1 parameter vector into named columns.

    Layout (see ``pFIC.get_parameters``): ``[wEE (N), wEI (N), G (1), sigma (N)]``.
    ``G`` is a single global coupling; it is repeated down the column so the file
    stays a plain rectangular table.
    """
    parameter = as_column(best_from_test["parameter"])
    expected = 3 * n_roi + 1
    if parameter.size != expected:
        raise SystemExit(
            f"[error] expected a {expected}-element parameter vector for "
            f"{n_roi} ROIs, got {parameter.size}.")
    frame = pd.DataFrame({
        "region_index": np.arange(n_roi),
        "wEE": parameter[0:n_roi],
        "wEI": parameter[n_roi:2 * n_roi],
        "sigma": parameter[2 * n_roi + 1:3 * n_roi + 1],
        "G": parameter[2 * n_roi],
    })
    path = os.path.join(out_dir, "best_parameters.csv")
    frame.to_csv(path, index=False)
    print(f"Wrote {path}")
    return frame


def write_fit_quality_csv(best_from_test, out_dir):
    """Report the FC+FCD cost of the chosen parameter set at every stage.

    ``ModelValidator`` carries each stage's losses forward with a stage prefix,
    so the final dict holds ``train_*``, ``val_*`` and un-prefixed test losses.
    """
    loss_names = ("corr_loss", "l1_loss", "ks_loss", "total_loss")
    rows = []
    for stage, prefix in (("train", "train_"), ("val", "val_"), ("test", "")):
        row = {"stage": stage}
        found = False
        for name in loss_names:
            key = f"{prefix}{name}"
            if key in best_from_test:
                row[name] = float(as_column(best_from_test[key])[0])
                found = True
            else:
                row[name] = np.nan
        if found:
            rows.append(row)
    frame = pd.DataFrame(rows)
    path = os.path.join(out_dir, "fit_quality.csv")
    frame.to_csv(path, index=False)
    print(f"Wrote {path}")
    return frame


def write_run_config(config, out_dir):
    path = os.path.join(out_dir, "run_config.ini")
    with open(path, "w") as handle:
        config.write(handle)
    print(f"Wrote {path}")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def estimate_ei(args):
    """Run the three CMA-ES stages for one subject and write the CSVs."""
    n_roi, n_bins = validate_inputs(args)
    if args.delssome_ckpt is not None:
        validate_ckpt(args.delssome_ckpt, n_roi, n_bins)

    set_torch_default()

    config = build_config(tr=args.tr,
                          scan_length=args.scan_length,
                          n_roi=n_roi,
                          base_ini=args.base_config,
                          overrides=args.config,
                          atlas=args.atlas)
    fcd_bins = int(config["Simulating Parameters"]["fcd_hist_bins"])
    if fcd_bins != n_bins:
        raise SystemExit(
            f"[error] the config expects {fcd_bins} FCD-CDF bins "
            f"(fcd_hist_bins) but --fcd-cdf has {n_bins} rows. Recompute your "
            "FCD CDF with matching bins, or override fcd_hist_bins via "
            "--config.")

    out_dir = os.path.abspath(args.out_dir)
    work_dir = os.path.join(out_dir, ".work")
    if args.overwrite and os.path.isdir(work_dir):
        shutil.rmtree(work_dir)
    stage_dirs = {stage: os.path.join(work_dir, stage) for stage in STAGES}
    for path in stage_dirs.values():
        os.makedirs(path, exist_ok=True)
    write_run_config(config, out_dir)

    # One emp_stats for all three stages. The DL-vectorised fields are only read
    # by the training stage, so carrying them through is harmless and saves
    # re-reading the CSVs twice.
    emp_stats = build_stats_from_paths(
        fc_path=args.fc,
        fcd_cdf_path=args.fcd_cdf,
        sc_path=args.sc,
        myelin_path=args.myelin,
        gradient_path=args.gradient,
        ds_name=args.subject_id,
        group_or_sub_id=args.subject_id,
        atlas=args.atlas,
        with_myelin_gradient=True,
        with_dl_stats=True,
    )
    use_delssome = args.delssome_ckpt is not None
    if use_delssome:
        emp_stats["delssome_ckpt"] = args.delssome_ckpt

    # --- stage 1: CMA-ES search ------------------------------------------- #
    print(f"\n=== stage 1/3: CMA-ES search "
          f"({'DELSSOME' if use_delssome else 'Euler'}, "
          f"{args.num_epochs} epochs) ===", flush=True)
    epochs = np.arange(0, args.num_epochs)
    state = train_helper(
        config=config,
        emp_stats=emp_stats,
        train_save_dir=stage_dirs["train"],
        num_epochs=args.num_epochs,
        dl_epoch_range=epochs if use_delssome else [],
        euler_epoch_range=[] if use_delssome else epochs,
        opportunities=args.opportunities,
        next_epoch=0,
        seed=args.seed,
    )
    if state != 0:
        print("[error] CMA-ES failed to produce a usable parameter population "
              f"(exit state {state}). With state 1 the search could not even "
              "initialise: check that your FC and FCD CDF are on the expected "
              "scales, and that TR and scan length are right.")
        return state
    get_best_train_params(stage_dirs["train"],
                          num_epochs=args.num_epochs,
                          top_k_for_each_epoch=1)

    # --- stage 2: validation (always Euler) -------------------------------- #
    print("\n=== stage 2/3: validation (Euler) ===", flush=True)
    ModelValidator(
        config, emp_stats,
        get_best_params_file_path("train", stage_dirs["train"]),
        stage_dirs["val"],
    ).validate(seed=args.seed)

    # --- stage 3: test (always Euler, finer dt) ---------------------------- #
    print("\n=== stage 3/3: test (Euler, finer dt) ===", flush=True)
    ModelTester(
        config, emp_stats,
        get_best_params_file_path("val", stage_dirs["val"]),
        stage_dirs["test"],
    ).test(seed=args.seed)

    # --- write the human-readable results --------------------------------- #
    best_from_test = torch.load(
        get_best_params_file_path("test", stage_dirs["test"]),
        map_location="cpu", weights_only=False,
    )
    print()
    ei = write_ei_csv(best_from_test, out_dir)
    write_parameters_csv(best_from_test, n_roi, out_dir)
    write_fit_quality_csv(best_from_test, out_dir)
    print(f"\nMean E/I ratio across {len(ei)} regions: "
          f"{ei['EI_ratio'].mean():.6g}")
    return 0


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Estimate a per-region E/I ratio for one subject from your "
        "own FC and FCD-CDF CSVs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = p.add_argument_group("empirical data (required)")
    required.add_argument("--fc", required=True,
                          help="Header-less N x N functional-connectivity CSV.")
    required.add_argument("--fcd-cdf", required=True,
                          help="Header-less bins x 1 FCD-CDF CSV.")
    required.add_argument("--tr", required=True, type=float,
                          help="Repetition time of your fMRI, in SECONDS. No "
                          "default on purpose: inheriting HCP-YA's 0.72 s "
                          "silently produces wrong results.")
    required.add_argument("--scan-length", required=True, type=float,
                          help="Length of your scan, in SECONDS "
                          "(e.g. 300 frames at TR 2.0 -> 600).")
    required.add_argument("--out-dir", required=True,
                          help="Directory for the result CSVs.")

    optional = p.add_argument_group("optional empirical data")
    optional.add_argument("--sc", default=None,
                          help="Header-less N x N structural-connectivity CSV. "
                          "Defaults to the bundled HCP-YA group SC, which is "
                          "what the paper uses for every subject.")
    optional.add_argument("--myelin", default=None,
                          help="Header-less N x 1 myelin CSV. Defaults to the "
                          "bundled HCP-YA group map. Must be given together "
                          "with --gradient.")
    optional.add_argument("--gradient", default=None,
                          help="Header-less N x 1 RSFC-gradient CSV. Defaults "
                          "to the bundled HCP-YA group map.")

    model = p.add_argument_group("optimisation")
    model.add_argument("--delssome-ckpt", default=None,
                       help="Trained DELSSOME predictor. When given, the search "
                       "stage uses it instead of Euler integration; validation "
                       "and test stay on Euler.")
    model.add_argument("--num-epochs", default=100, type=int,
                       help="CMA-ES epochs (100 candidates each).")
    model.add_argument("--seed", default=None, type=int,
                       help="RNG seed. Set it for reproducible output.")
    model.add_argument("--opportunities", default=10, type=int,
                       help="How many times to restart CMA-ES before giving up.")
    model.add_argument("--atlas", default="DK68",
                       help="Only used to locate the bundled default SC/myelin/"
                       "gradient maps.")
    model.add_argument("--subject-id", default="subject",
                       help="Label used in log messages only.")
    model.add_argument("--base-config", default=None,
                       help="INI supplying the biophysical and hemodynamic "
                       "constants. Defaults to the released "
                       "configs/model/pfic/<atlas>/indiv/config_HCP-YA.ini; "
                       "those constants are dataset-independent.")
    model.add_argument("--config", default=None,
                       help="Optional second INI layered on top, to change e.g. "
                       "dt_train or param_dup.")
    model.add_argument("--overwrite", action="store_true",
                       help="Delete <out-dir>/.work before running.")
    return p


def main():
    args = build_parser().parse_args()
    print(datetime.datetime.now(), ": MAIN PROGRAM START", flush=True)
    print(f"Default configs are read from: {CONFIG_DIR}")
    if torch.cuda.is_available():
        print("Current GPU:",
              torch.cuda.get_device_name(torch.cuda.current_device()))
    state = estimate_ei(args)
    print(datetime.datetime.now(), ": MAIN PROGRAM END", flush=True)
    sys.exit(state)


if __name__ == "__main__":
    main()
