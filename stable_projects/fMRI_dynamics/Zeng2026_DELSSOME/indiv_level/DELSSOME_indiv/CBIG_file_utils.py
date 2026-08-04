# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Path builders, empirical-data loaders, and parameter-dict combiners.

Only the functionality required by the individual-level DELSSOME pipeline is
included here: the per-subject (group-type ``None``) code paths. ComBat and
group-level helpers are not part of this release.

The expected directory structure is:

    LOG_DIR/<ds>/<target>/<phase>/trial<N>/seed<N>/<sub_id>/
        param_save_epoch<k>.pth   (one per training epoch)
        best_from_train.pth
        best_from_val.pth         (after validation)
        ...

    DATA_DIR/<ds>/FC/<atlas>/<sub_id>_bld00{0..6}.csv
    DATA_DIR/<ds>/FCD_cdf/<atlas>/<sub_id>_bld00{0..6}.csv
    DATA_DIR/HCP-YA/pFIC_input/<atlas>/{SC_train_group_dl_ds, myelin, rsfc_gradient}.csv
    DATA_DIR/<ds>/demogr/demogr.csv

Override ``DELSSOME_DATA_DIR`` / ``DELSSOME_LOG_DIR`` / ``DELSSOME_CONFIG_DIR``
environment variables to point at your own data location.
"""

import configparser
import os
import numpy as np
import pandas as pd
import torch

from DELSSOME_indiv.CBIG_constants import (
    CONFIG_DIR,
    DATA_DIR,
    LOG_DIR,
    DEFAULT_DTYPE,
    HCPYA_pFIC_INPUT,
    PREV_PHASE,
)
from DELSSOME_indiv.CBIG_FC_FCD_utils import fisher_average


# ========================================================================== #
# Path builders
# ========================================================================== #
def get_target_dir(ds_name, target):
    path = os.path.join(LOG_DIR, ds_name, target)
    os.makedirs(path, exist_ok=True)
    return path


def get_phase_dir(ds_name, target, phase):
    path = os.path.join(get_target_dir(ds_name, target), phase)
    os.makedirs(path, exist_ok=True)
    return path


def get_trial_dir(ds_name, target, phase, trial_idx):
    path = os.path.join(get_phase_dir(ds_name, target, phase), f"trial{trial_idx}")
    os.makedirs(path, exist_ok=True)
    return path


def get_run_dir(ds_name, target, phase, trial_idx, seed_idx, group_or_sub_id=None):
    """Returns ``<LOG_DIR>/<ds>/<target>/<phase>/trial<N>/seed<N>[/<sub_id>]``."""
    run_path = os.path.join(
        get_target_dir(ds_name, target),
        phase,
        f"trial{trial_idx}",
        f"seed{seed_idx}",
    )
    if group_or_sub_id is not None:
        run_path = os.path.join(run_path, str(group_or_sub_id))
    os.makedirs(run_path, exist_ok=True)
    return run_path


# ========================================================================== #
# Per-epoch / per-phase file names
# ========================================================================== #
def get_train_epoch_file_path(train_save_dir, epoch_idx):
    return os.path.join(train_save_dir, f"param_save_epoch{epoch_idx}.pth")


def get_best_params_file_path(phase, save_dir):
    """Path of the bootstrap params for the NEXT phase (e.g. ``best_from_train.pth``)."""
    return os.path.join(save_dir, f"best_from_{phase}.pth")


def get_agg_seeds_range(agg_seeds_num, agg_seed_idx):
    """Map an aggregated-seed index to its underlying seed range + label.

    Useful when validation/test merge ``agg_seeds_num`` consecutive training
    seeds into one parameter pool (paper Section 2.3).
    """
    agg_seeds_num = int(agg_seeds_num)
    agg_seed_idx = int(agg_seed_idx)
    seed_start = (agg_seed_idx - 1) * agg_seeds_num + 1
    seed_end = agg_seed_idx * agg_seeds_num
    agg_seed_label = f"_seed{seed_start}_to_seed{seed_end}"
    return range(seed_start, seed_end + 1), agg_seed_label


def get_curr_phase_save_dir(ds_name, target, phase, trial_idx, seed_idx,
                            group_or_sub_id=None, agg_seeds_num=None):
    if agg_seeds_num is not None and phase != "train":
        _, seed_idx = get_agg_seeds_range(agg_seeds_num, seed_idx)
    return get_run_dir(ds_name, target, phase, trial_idx, seed_idx,
                       group_or_sub_id=group_or_sub_id)


def get_prev_phase_best_params_path(ds_name, target, curr_phase, trial_idx,
                                    seed_idx, curr_group_or_sub_id=None,
                                    group_type=None, agg_seeds_num=None):
    """Locate the best params dumped by the phase BEFORE ``curr_phase``."""
    if agg_seeds_num is not None and curr_phase != "train":
        seeds_range, seed_idx = get_agg_seeds_range(agg_seeds_num, seed_idx)
    prev_phase = PREV_PHASE[curr_phase]
    prev_group_or_sub_id = curr_group_or_sub_id
    if group_type is not None:
        prev_group_or_sub_id = curr_group_or_sub_id.replace(curr_phase, prev_phase)

    prev_phase_save_dir = get_run_dir(
        ds_name, target, prev_phase, trial_idx, seed_idx,
        group_or_sub_id=prev_group_or_sub_id,
    )
    prev_phase_best_params_path = get_best_params_file_path(
        prev_phase, prev_phase_save_dir
    )

    # If best params from individual seeds were not yet aggregated, do it now.
    if (agg_seeds_num is not None
            and curr_phase != "train"
            and not os.path.exists(prev_phase_best_params_path)):
        prev_phase_all_seeds_save_dir = [
            get_run_dir(ds_name, target, prev_phase, trial_idx, seed_i,
                        prev_group_or_sub_id)
            for seed_i in seeds_range
        ]
        prev_phase_all_seeds_best_params_path = [
            get_best_params_file_path(prev_phase, d)
            for d in prev_phase_all_seeds_save_dir
        ]
        combine_all_param_dicts(
            paths_to_dicts=prev_phase_all_seeds_best_params_path,
            combined_dict_save_path=prev_phase_best_params_path,
        )
    return prev_phase_best_params_path


# ========================================================================== #
# Parameter-dict combiners (used to build best_from_train.pth)
# ========================================================================== #
def get_values_at_indices(saved_dict: dict, indices=None):
    """Take ``saved_dict[key][indices]`` (or the last-axis equivalent for 2D)."""
    if indices is None:
        return saved_dict
    new_dict = {}
    for key, value in saved_dict.items():
        if not torch.is_tensor(value):
            value = torch.as_tensor(value)
        if value.ndim == 0:
            new_dict[key] = value
            continue
        assert value.ndim in (1, 2), f"The value for {key} is not a 1D or 2D tensor"
        value = value[indices] if value.ndim == 1 else value[:, indices]
        new_dict[key] = value
    return new_dict


def sort_dict_by_total_loss(saved_dict, top_k=None):
    sorted_indices = torch.argsort(saved_dict["total_loss"], descending=False)
    new_save_dict = get_values_at_indices(saved_dict, sorted_indices)
    if top_k is not None:
        new_save_dict = get_values_at_indices(new_save_dict, range(top_k))
    return new_save_dict


def get_first_k_values(values, k=None):
    """Take the first ``k`` columns of a 2D tensor (or values of a 1D tensor)."""
    if not torch.is_tensor(values):
        values = torch.as_tensor(values)
    if values.ndim == 0:
        return values.unsqueeze(0).unsqueeze(0)
    if values.ndim == 1:
        if k is not None:
            values = values[:k]
        values = values.unsqueeze(0)
    elif values.ndim == 2:
        if k is not None:
            values = values[:, :k]
    else:
        raise Exception("The input `values` is not a 1D or 2D tensor")
    return values


def combine_all_param_dicts(paths_to_dicts,
                            are_params_sorted=True,
                            top_k_per_dict=None,
                            top_k_among_all_dicts=None,
                            combined_dict_save_path=None,
                            ignore_keys=None,
                            dict_ds_indices=None,
                            dict_group_or_sub_indices=None):
    """
    Merge a list of per-epoch (or per-seed) param dictionaries into one.

    Keeps only the top-k parameter sets per dictionary, optionally adds dataset
    and subject indices, and finally sorts the merged dictionary by total loss.
    """
    ignore_keys = ignore_keys or []
    combined_dict = {}
    for i, path_to_dict in enumerate(paths_to_dicts):
        if not os.path.exists(path_to_dict):
            print(f"[Warning]: File {path_to_dict} does not exist!")
            continue
        d = torch.load(path_to_dict, map_location="cpu")
        if not are_params_sorted:
            d["parameter"] = d["parameter"][:, d["valid_param_indices"]]
        if dict_ds_indices is not None:
            d["ds_indices"] = torch.full(d["parameter"].shape[1:],
                                         dict_ds_indices[i])
        if dict_group_or_sub_indices is not None:
            d["group_or_sub_indices"] = torch.full(
                d["parameter"].shape[1:], dict_group_or_sub_indices[i])
        for key, value in d.items():
            if key in ignore_keys:
                continue
            combined_dict.setdefault(key, [])
            value = get_first_k_values(value, k=top_k_per_dict)
            combined_dict[key].append(value)

    if not combined_dict:
        raise Exception("All input param dict paths do not exist!")

    for key, value in combined_dict.items():
        combined_dict[key] = torch.cat(value, dim=1).squeeze()

    combined_dict = sort_dict_by_total_loss(combined_dict, top_k_among_all_dicts)

    if combined_dict_save_path is not None:
        os.makedirs(os.path.dirname(combined_dict_save_path), exist_ok=True)
        torch.save(combined_dict, combined_dict_save_path)
        print(f"Combined dict saved to {combined_dict_save_path}")
    return combined_dict


def get_best_train_params(train_save_dir, num_epochs=100, top_k_for_each_epoch=1):
    """
    Build ``best_from_train.pth`` by aggregating the top-k candidates from
    every training epoch (one row per epoch by default).
    """
    train_save_files = [
        get_train_epoch_file_path(train_save_dir, ep) for ep in range(num_epochs)
    ]
    best_from_train_file_path = get_best_params_file_path("train", train_save_dir)
    if top_k_for_each_epoch != 1:
        best_from_train_file_path = best_from_train_file_path.replace(
            ".pth", f"_top{top_k_for_each_epoch}_per_epoch.pth"
        )
    best_from_train = combine_all_param_dicts(
        are_params_sorted=False,
        paths_to_dicts=train_save_files,
        top_k_per_dict=top_k_for_each_epoch,
        combined_dict_save_path=best_from_train_file_path,
    )
    print(
        f"Successfully saved the top {top_k_for_each_epoch} parameters from each "
        f"train epoch to: {best_from_train_file_path}"
    )
    return best_from_train


# ========================================================================== #
# Empirical-data loaders
# ========================================================================== #
def get_demogr(ds_name, group_type=None):
    """Load the per-dataset demographics table (``demogr.csv``)."""
    if group_type is None:
        demogr_file_path = os.path.join(DATA_DIR, ds_name, "demogr", "demogr.csv")
    else:
        demogr_file_path = os.path.join(DATA_DIR, ds_name, "demogr", group_type,
                                        "demogr.csv")
    if not os.path.exists(demogr_file_path):
        raise FileNotFoundError(f"Missing demographics file: {demogr_file_path}")
    return pd.read_csv(demogr_file_path)


def load_matrix_csv(path, dtype=DEFAULT_DTYPE, what="matrix"):
    """Read a header-less numeric CSV into a 2D tensor, with a useful error."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {what}: {path}")
    values = pd.read_csv(path, header=None).values.copy()
    return torch.as_tensor(values).to(dtype)


def get_myelin(ds_name, group_or_sub_id, group_type, atlas="DK68",
               dtype=DEFAULT_DTYPE, path=None):
    """Load the regional myelin profile.

    All datasets in the paper share the HCP-YA group myelin map, so ``ds_name``
    and ``group_or_sub_id`` are accepted for signature symmetry but unused when
    ``path`` is None. Pass ``path`` to supply your own N-by-1 CSV.
    """
    emp_myelin_path = path or os.path.join(HCPYA_pFIC_INPUT, atlas, "myelin.csv")
    print("Loading myelin from:", emp_myelin_path)
    return load_matrix_csv(emp_myelin_path, dtype=dtype, what="myelin CSV")


def get_rsfc_gradient(ds_name, group_or_sub_id, group_type, atlas="DK68",
                      dtype=DEFAULT_DTYPE, path=None):
    """Load the regional RSFC-gradient profile.

    Defaults to the HCP-YA group map shared by every dataset in the paper; pass
    ``path`` to supply your own N-by-1 CSV.
    """
    emp_path = path or os.path.join(HCPYA_pFIC_INPUT, atlas, "rsfc_gradient.csv")
    print("Loading rsfc gradient from:", emp_path)
    return load_matrix_csv(emp_path, dtype=dtype, what="RSFC-gradient CSV")


def get_sc_mat(ds_name, group_or_sub_id, group_type=None, atlas="DK68",
               dtype=DEFAULT_DTYPE, path=None):
    """Load the structural-connectivity matrix.

    The paper uses the HCP-YA training-set group SC matrix for every subject,
    since per-subject SC was not available across all 14 datasets. To use your
    own (group or per-subject) SC, pass ``path`` — no source edit needed.
    """
    emp_sc_mat_path = path or os.path.join(HCPYA_pFIC_INPUT, atlas,
                                           "SC_train_group_dl_ds.csv")
    print("Loading SC from:", emp_sc_mat_path)
    return load_matrix_csv(emp_sc_mat_path, dtype=dtype, what="SC CSV")


def get_emp_fc(ds_name, group_or_sub_id, group_type, atlas="DK68",
               use_harmonized=False, dtype=DEFAULT_DTYPE, verbose=True,
               path=None):
    """Load the empirical FC matrix for a subject (Fisher-averaged over runs).

    Pass ``path`` to read a single explicit CSV instead of resolving
    ``DATA_DIR/<ds>/FC/<atlas>/<sub>_bld00{0..6}.csv`` by convention.
    """
    if path is not None:
        emp_fc_paths = [path]
    else:
        fc_dir = "harmonized_FC" if use_harmonized else "FC"
        average_fc_across_runs_path = os.path.join(
            DATA_DIR, ds_name, fc_dir, atlas, f"{group_or_sub_id}_bld000.csv"
        )
        if os.path.exists(average_fc_across_runs_path):
            emp_fc_paths = [average_fc_across_runs_path]
        else:
            emp_fc_paths = [
                os.path.join(DATA_DIR, ds_name, fc_dir, atlas,
                             f"{group_or_sub_id}_bld00{i}.csv")
                for i in range(1, 7)
            ]
    if verbose:
        print("Loading FC from:", emp_fc_paths)
    emp_fcs = [pd.read_csv(p, header=None).values.copy()
               for p in emp_fc_paths if os.path.exists(p)]
    if not emp_fcs:
        raise ValueError(
            f"Empirical FC data not found for subject/group {group_or_sub_id}. "
            f"Looked in: {emp_fc_paths}"
        )
    emp_fc = fisher_average(np.array(emp_fcs)) if len(emp_fcs) > 1 else emp_fcs[0]
    return torch.as_tensor(emp_fc).to(dtype)


def get_emp_fcd_cdf(ds_name, group_or_sub_id, group_type, atlas="DK68",
                    use_harmonized=False, dtype=DEFAULT_DTYPE,
                    verbose=True, normalize=True, path=None):
    """Load the empirical FCD CDF (averaged across runs, normalised to peak 1).

    Pass ``path`` to read a single explicit CSV instead of resolving
    ``DATA_DIR/<ds>/FCD_cdf/<atlas>/<sub>_bld00{0..6}.csv`` by convention.
    """
    if path is not None:
        emp_fcd_cdf_paths = [path]
    else:
        fcd_dir = "harmonized_FCD_cdf" if use_harmonized else "FCD_cdf"
        avg_path = os.path.join(
            DATA_DIR, ds_name, fcd_dir, atlas, f"{group_or_sub_id}_bld000.csv"
        )
        if os.path.exists(avg_path):
            emp_fcd_cdf_paths = [avg_path]
        else:
            emp_fcd_cdf_paths = [
                os.path.join(DATA_DIR, ds_name, fcd_dir, atlas,
                             f"{group_or_sub_id}_bld00{i}.csv")
                for i in range(1, 7)
            ]
    if verbose:
        print("Loading FCD CDF from:", emp_fcd_cdf_paths)
    emp_fcds = [pd.read_csv(p, header=None).values.flatten()
                for p in emp_fcd_cdf_paths if os.path.exists(p)]
    if not emp_fcds:
        raise ValueError(
            f"Empirical FCD data not found for subject/group {group_or_sub_id}. "
            f"Looked in: {emp_fcd_cdf_paths}"
        )
    emp_fcd = np.mean(emp_fcds, axis=0) if len(emp_fcds) > 1 else emp_fcds[0]
    emp_fcd = torch.as_tensor(emp_fcd).to(dtype)
    if normalize:
        emp_fcd = emp_fcd / emp_fcd[-1]
    emp_fcd = emp_fcd.unsqueeze(1)
    assert emp_fcd.ndim == 2 and emp_fcd.shape[1] == 1, (
        f"Empirical FCD CDF should be a column vector, got shape {emp_fcd.shape}"
    )
    return emp_fcd


def get_dl_stats(stats):
    """Pre-compute the DL-friendly vectorised versions of SC/FC and the FCD pdf.

    Used by the DELSSOME cost predictor at training time and inside DELSSOME
    CMA-ES (apply phase) at inference time.
    """
    N = stats["sc_mat"].shape[0]
    upper_trig_mask = torch.triu(torch.ones(N, N, dtype=torch.bool), 1)
    sc_dl = stats["sc_mat"][upper_trig_mask]   # [N*(N-1)/2]
    fc_dl = stats["emp_fc"][upper_trig_mask]   # [N*(N-1)/2]
    # FCD-PDF is approximated from the CDF: scale by 100 (paper convention).
    fcd_dl = torch.diff(
        stats["emp_fcd_cdf"].squeeze() * 100,
        dim=0,
        prepend=torch.as_tensor([0]),
    )
    assert fcd_dl.ndim == 1
    stats["sc_dl"] = sc_dl
    stats["fc_dl"] = fc_dl
    stats["fcd_dl"] = fcd_dl
    return stats


def build_stats_from_paths(fc_path, fcd_cdf_path, sc_path=None,
                           myelin_path=None, gradient_path=None,
                           ds_name="user", group_or_sub_id="subject",
                           group_type=None, atlas="DK68",
                           with_myelin_gradient=True, with_dl_stats=False,
                           use_harmonized=False, dtype=DEFAULT_DTYPE,
                           verbose=True):
    """
    Build the ``emp_stats`` dict from explicit file paths.

    This is the single place where empirical inputs become the dict that
    ``CmaesTrainer`` / ``ModelValidator`` / ``ModelTester`` and ``pFIC`` consume:

      * ``sc_mat`` — raw structural connectivity matrix.
      * ``sc_euler`` — SC scaled by ``0.02 / max(SC)`` (used inside Euler integration).
      * ``emp_fc``, ``emp_fcd_cdf`` — empirical FC and FCD CDF.
      * Optional: ``myelin``, ``rsfc_gradient`` (parameterisation regressors).
      * Optional: ``sc_dl``, ``fc_dl``, ``fcd_dl`` (DL-friendly vectorised form).

    ``sc_path`` / ``myelin_path`` / ``gradient_path`` may be None, in which case
    the bundled HCP-YA group maps under ``HCPYA_pFIC_INPUT`` are used — the same
    behaviour as the paper. :func:`get_stats` resolves paths by the on-disk
    naming convention and then delegates here, so both entry points produce
    byte-identical dicts.
    """
    sc_mat = get_sc_mat(ds_name=ds_name, group_or_sub_id=group_or_sub_id,
                        group_type=group_type, atlas=atlas, dtype=dtype,
                        path=sc_path)
    emp_fc = get_emp_fc(ds_name=ds_name, group_or_sub_id=group_or_sub_id,
                        group_type=group_type, atlas=atlas,
                        use_harmonized=use_harmonized, dtype=dtype,
                        verbose=verbose, path=fc_path)
    # Normalised SC used by Euler integration of the FIC ODEs.
    sc_euler = sc_mat / torch.max(sc_mat) * 0.02
    emp_fcd_cdf = get_emp_fcd_cdf(
        ds_name=ds_name, group_or_sub_id=group_or_sub_id,
        group_type=group_type, atlas=atlas,
        use_harmonized=use_harmonized, dtype=dtype,
        verbose=verbose, path=fcd_cdf_path,
    )
    stats = {
        "ds_name": ds_name,
        "group_or_sub_id": group_or_sub_id,
        "group_type": group_type,
        "atlas": atlas,
        "sc_mat": sc_mat,
        "sc_euler": sc_euler,
        "emp_fc": emp_fc,
        "emp_fcd_cdf": emp_fcd_cdf,
    }
    if with_myelin_gradient:
        stats["myelin"] = get_myelin(ds_name, group_or_sub_id, group_type,
                                     atlas=atlas, dtype=dtype,
                                     path=myelin_path)
        stats["rsfc_gradient"] = get_rsfc_gradient(
            ds_name, group_or_sub_id, group_type, atlas=atlas, dtype=dtype,
            path=gradient_path,
        )
    if with_dl_stats:
        stats = get_dl_stats(stats)
    print("Statistics for subject/group", group_or_sub_id, "loaded")
    return stats


def get_stats(ds_name, group_or_sub_id, group_type=None, atlas="DK68",
              with_myelin_gradient=True, with_dl_stats=False,
              use_harmonized=False, dtype=DEFAULT_DTYPE):
    """
    Collect every empirical input the FIC model needs for one subject, resolving
    paths from the on-disk convention documented at the top of this module.

    Thin wrapper over :func:`build_stats_from_paths` with ``*_path=None``, i.e.
    "resolve everything by convention". Use ``build_stats_from_paths`` directly
    when the user supplies explicit CSVs.
    """
    return build_stats_from_paths(
        fc_path=None, fcd_cdf_path=None, sc_path=None,
        myelin_path=None, gradient_path=None,
        ds_name=ds_name, group_or_sub_id=group_or_sub_id,
        group_type=group_type, atlas=atlas,
        with_myelin_gradient=with_myelin_gradient,
        with_dl_stats=with_dl_stats, use_harmonized=use_harmonized,
        dtype=dtype,
    )


# ========================================================================== #
# INI configuration
# ========================================================================== #
def read_config_file(config_path):
    """Read an INI config, failing loudly when the file is absent.

    ``ConfigParser.read`` silently ignores a missing path, which used to turn a
    typo in ``--ds-name`` into a confusing ``KeyError: 'Dynamic Model Constants'``
    thousands of lines later. Fail here instead, naming the path we looked for.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"Config file not found: {config_path}\n"
            "Individual-level configs are resolved as\n"
            "  CONFIG_DIR/model/<dynamic_model>/<atlas>/indiv/config_<ds_name>.ini\n"
            "so a new dataset needs its own config there. Copy\n"
            "  configs/model/pfic/DK68/indiv/config_template.ini\n"
            "and set TR / simulate_time / window_size for your scans."
        )
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def apply_timing_to_config(config, tr, scan_length):
    """Set the four scan-timing keys from a TR and a scan length (seconds).

    ``window_size`` and ``burn_in_time`` are derived exactly as the paper does:
    a ~60 s FCD window and a 72 s burn-in, both rounded to whole TRs.
    """
    tr = float(tr)
    scan_length = float(scan_length)
    if tr <= 0:
        raise ValueError(f"TR must be positive, got {tr}")
    if scan_length <= 0:
        raise ValueError(f"scan_length must be positive, got {scan_length}")
    config["Dataset Parameters"]["TR"] = str(tr)
    config["Dataset Parameters"]["simulate_time"] = str(scan_length)
    config["Dataset Parameters"]["window_size"] = str(round(60 / tr))
    config["Dataset Parameters"]["burn_in_time"] = str(np.ceil(72 / tr) * tr)
    return config


def build_config(tr, scan_length, n_roi=None, base_ini=None, overrides=None,
                 atlas="DK68", dynamic_model_name="pfic"):
    """
    Build a config from an explicit TR / scan length, with no dataset lookup.

    This is the path used by ``CBIG_DELSSOME_indiv_estimate_EI.py``: the user
    states their own scan timing rather than silently inheriting HCP-YA's
    0.72 s TR / 864 s scan length from ``config_HCP-YA.ini``.

    Args:
        tr: repetition time of the empirical fMRI, in seconds.
        scan_length: length of the empirical scan, in seconds.
        n_roi: number of ROIs; when given, overrides the base config.
        base_ini: INI supplying the biophysical / hemodynamic constants.
            Defaults to the released HCP-YA config, whose ``[Dynamic Model
            Constants]`` and ``[Hemodynamic Model Constants]`` are atlas- and
            dataset-independent.
        overrides: optional path to a second INI whose keys are layered on top,
            for changing e.g. ``dt_train`` or ``param_dup``.
    """
    if base_ini is None:
        base_ini = os.path.join(CONFIG_DIR, "model", dynamic_model_name, atlas,
                                "indiv", "config_HCP-YA.ini")
    config = read_config_file(base_ini)
    config = apply_timing_to_config(config, tr, scan_length)
    if n_roi is not None:
        config["Simulating Parameters"]["n_ROI"] = str(int(n_roi))
    if overrides is not None:
        extra = read_config_file(overrides)
        for section in extra.sections():
            if not config.has_section(section):
                config.add_section(section)
            for key, value in extra.items(section):
                config[section][key] = value
    return config


def get_config(ds_name, atlas, group_or_sub_id=None, group_type=None,
               dynamic_model_name="pfic"):
    """
    Load the per-subject INI config (Dataset / Simulating / Dynamic / Hemodynamic).

    The base config lives at::

        CONFIG_DIR/model/<dynamic_model_name>/<atlas>/<indiv|group_type>/config_<ds>.ini

    If the dataset ships a ``demogr.csv`` carrying ``TR`` and ``scan_length``,
    those override the INI's scan-timing keys -- for any dataset, not just the
    ones used in the paper. Datasets whose ``demogr.csv`` lacks those columns
    keep the INI values, so check them before running your own cohort: the
    wrong scan timing produces a wrong E/I ratio rather than an error.
    """
    group_or_indiv_folder = "indiv" if group_type is None else group_type
    config_path = os.path.join(
        CONFIG_DIR, "model", dynamic_model_name, atlas,
        group_or_indiv_folder, f"config_{ds_name}.ini",
    )
    config = read_config_file(config_path)

    if group_or_sub_id is None:
        return config

    id_column = "group_or_sub_id" if group_type is not None else "subject_id"
    try:
        demogr = get_demogr(ds_name=ds_name, group_type=group_type)
    except FileNotFoundError:
        # No demographics table: the INI's own scan timing is authoritative.
        return config
    if not {"TR", "scan_length", id_column}.issubset(demogr.columns):
        return config
    row = demogr[demogr[id_column] == group_or_sub_id]
    if row.empty:
        print(f"[warn] {group_or_sub_id} not in {ds_name} demogr.csv; "
              f"using the scan timing from {os.path.basename(config_path)}")
        return config
    return apply_timing_to_config(config, row["TR"].values[0],
                                  row["scan_length"].values[0])


# ========================================================================== #
# Misc
# ========================================================================== #
def load_best_params(ds_name, target, phase, trial_idx, seed_idx,
                     group_or_sub_id=None, agg_seeds_num=None):
    """Load ``best_from_<phase>.pth`` if it exists, else ``None``."""
    curr_phase_save_dir = get_curr_phase_save_dir(
        ds_name, target, phase, trial_idx, seed_idx,
        group_or_sub_id, agg_seeds_num,
    )
    best_params_file_path = get_best_params_file_path(phase, curr_phase_save_dir)
    if not os.path.exists(best_params_file_path):
        return None
    return torch.load(best_params_file_path, map_location="cpu")
