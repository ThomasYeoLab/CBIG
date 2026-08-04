# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Collect per-subject E/I ratios into the one table the GAMLSS step wants.

``CBIG_DELSSOME_indiv_estimate_EI.py`` writes one ``EI_ratio.csv`` per subject.
``GAMLSS/CBIG_fit_gamlss_model.R`` wants a single table with one ROW per subject
and the columns ``age_in_years``, ``numeric_sex``, ``site_scanner`` and
``ei_ratio_mean``. This script bridges the two.

Example::

    cd indiv_level
    export PYTHONPATH=$PWD

    # after running estimate_EI once per subject into results/<sub_id>/
    python -m scripts.CBIG_DELSSOME_indiv_collect_EI \\
        --in-dir  results \\
        --demogr  my_data/demogr.csv \\
        --out     gamlss_input.csv

It also accepts the replication layout: pass ``--from-pth`` to read
``best_from_test.pth`` files (as written by ``CBIG_DELSSOME_indiv_cmaes.py``)
instead of ``EI_ratio.csv``.

The demographics table needs a subject-ID column plus whatever the GAMLSS
formula uses. With the defaults that means:

    subject_id,age_in_years,sex,site,scanner
    sub-000,34.0,female,WashU,Siemens

``numeric_sex`` is derived as ``sex == 'female'`` (0 = male, 1 = female) when not
already present, and ``site_scanner`` as ``"<site>_<scanner>"``. Both are used
verbatim if your table already provides them.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import torch

#: Column carrying the mean cortical E/I ratio -- the GAMLSS default response.
MEAN_COLUMN = "ei_ratio_mean"
#: Per-region columns, so a per-ROI GAMLSS can be fitted via --target.
ROI_COLUMN_FMT = "ei_ratio_roi_{}"


def find_subject_dirs(in_dir, filename):
    """Return ``[(subject_id, path), ...]`` for every subject under ``in_dir``.

    ``subject_id`` is the directory name holding the result file, which is how
    ``estimate_EI --out-dir results/<sub_id>`` and the replication layout
    (``.../<sub_id>/best_from_test.pth``) both name their subjects.
    """
    found = []
    for root, _dirs, files in os.walk(in_dir):
        if filename in files:
            found.append((os.path.basename(root), os.path.join(root, filename)))
    return sorted(found)


def read_ei_from_csv(path):
    frame = pd.read_csv(path)
    if "EI_ratio" not in frame.columns:
        raise SystemExit(f"[error] {path} has no 'EI_ratio' column.")
    return frame["EI_ratio"].to_numpy(dtype=float)


def read_ei_from_pth(path):
    saved = torch.load(path, map_location="cpu", weights_only=False)
    if "EI_for_valid_params" not in saved:
        raise SystemExit(
            f"[error] {path} has no 'EI_for_valid_params'. Only the 'test' "
            "stage writes E/I estimates.")
    values = torch.as_tensor(saved["EI_for_valid_params"]).detach().cpu().numpy()
    # Shape is [N_roi, n_valid_params]; the test stage keeps a single winner.
    return np.asarray(values).reshape(values.shape[0], -1)[:, 0].astype(float)


def derive_gamlss_columns(demogr):
    """Add ``numeric_sex`` and ``site_scanner`` if the table lacks them."""
    demogr = demogr.copy()
    if "numeric_sex" not in demogr.columns:
        if "sex" not in demogr.columns:
            raise SystemExit(
                "[error] demographics table needs either 'numeric_sex' or "
                "'sex'. GAMLSS models sex as a fixed effect.")
        sex = demogr["sex"].astype(str).str.strip().str.lower()
        known = sex.isin(["male", "female", "m", "f"])
        if not known.all():
            raise SystemExit(
                "[error] 'sex' must be male/female (or m/f); found "
                f"{sorted(set(sex[~known]))}. Add a 'numeric_sex' column "
                "(0 = male, 1 = female) if your coding differs.")
        demogr["numeric_sex"] = sex.isin(["female", "f"]).astype(int)
    if "site_scanner" not in demogr.columns:
        if {"site", "scanner"}.issubset(demogr.columns):
            demogr["site_scanner"] = (demogr["site"].astype(str) + "_"
                                      + demogr["scanner"].astype(str))
        else:
            # A single-site cohort still needs the column: GAMLSS fits
            # random(site_scanner), which degenerates harmlessly to one level.
            print("[warn] no 'site'/'scanner' columns; setting site_scanner to "
                  "'single_site'. With one level the GAMLSS random effect does "
                  "no harmonisation.")
            demogr["site_scanner"] = "single_site"
    return demogr


def collect(args):
    filename = "best_from_test.pth" if args.from_pth else "EI_ratio.csv"
    reader = read_ei_from_pth if args.from_pth else read_ei_from_csv
    subjects = find_subject_dirs(args.in_dir, filename)
    if not subjects:
        raise SystemExit(
            f"[error] no '{filename}' found anywhere under {args.in_dir}. "
            "Run CBIG_DELSSOME_indiv_estimate_EI.py first, one --out-dir per "
            "subject (or pass --from-pth for the replication layout).")
    print(f"Found {len(subjects)} subject(s) with {filename} under {args.in_dir}")

    rows = []
    n_roi = None
    for subject_id, path in subjects:
        ei = reader(path)
        if n_roi is None:
            n_roi = ei.size
        elif ei.size != n_roi:
            raise SystemExit(
                f"[error] {subject_id} has {ei.size} regions but earlier "
                f"subjects had {n_roi}. All subjects must share a parcellation.")
        row = {"subject_id": subject_id, MEAN_COLUMN: float(np.mean(ei))}
        if not args.no_roi_columns:
            for idx, value in enumerate(ei):
                row[ROI_COLUMN_FMT.format(idx)] = float(value)
        rows.append(row)
    table = pd.DataFrame(rows)

    if args.demogr is not None:
        if not os.path.exists(args.demogr):
            raise SystemExit(f"[error] --demogr: file not found: {args.demogr}")
        demogr = derive_gamlss_columns(pd.read_csv(args.demogr))
        if args.id_column not in demogr.columns:
            raise SystemExit(
                f"[error] --demogr has no '{args.id_column}' column "
                f"(found: {list(demogr.columns)}). Use --id-column to name it.")
        demogr = demogr.rename(columns={args.id_column: "subject_id"})
        demogr["subject_id"] = demogr["subject_id"].astype(str)
        table["subject_id"] = table["subject_id"].astype(str)
        merged = table.merge(demogr, on="subject_id", how="left")
        missing = merged[merged["age_in_years"].isna()]["subject_id"].tolist() \
            if "age_in_years" in merged.columns else []
        if missing:
            print(f"[warn] {len(missing)} subject(s) had no demographics row "
                  f"and will be dropped by GAMLSS: {missing[:5]}"
                  f"{' ...' if len(missing) > 5 else ''}")
        # Put the modelling columns first so the file is readable by eye.
        lead = [c for c in ("subject_id", "age_in_years", "numeric_sex",
                            "site_scanner", "sex", "site", "scanner",
                            MEAN_COLUMN) if c in merged.columns]
        table = merged[lead + [c for c in merged.columns if c not in lead]]
    else:
        print("[warn] no --demogr given: the output has E/I ratios only. "
              "GAMLSS additionally needs age_in_years, numeric_sex and "
              "site_scanner.")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    table.to_csv(args.out, index=False)
    print(f"Wrote {args.out}: {len(table)} subjects x {len(table.columns)} columns")
    if len(table) < 100:
        print(f"[note] {len(table)} subjects is far below the ~12,000 used in "
              "the paper. GAMLSS needs a few hundred at minimum for a "
              "meaningful normative trajectory.")
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        description="Collect per-subject E/I ratios into a GAMLSS input table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--in-dir", required=True,
                   help="Directory searched recursively for per-subject "
                   "results. The directory NAME becomes the subject ID.")
    p.add_argument("--out", required=True, help="Output CSV path.")
    p.add_argument("--demogr", default=None,
                   help="Demographics CSV to merge in (age / sex / site / "
                   "scanner). Strongly recommended: GAMLSS needs it.")
    p.add_argument("--id-column", default="subject_id",
                   help="Subject-ID column in the demographics CSV.")
    p.add_argument("--from-pth", action="store_true",
                   help="Read best_from_test.pth instead of EI_ratio.csv, for "
                   "output produced by CBIG_DELSSOME_indiv_cmaes.py.")
    p.add_argument("--no-roi-columns", action="store_true",
                   help="Write only ei_ratio_mean, omitting the per-ROI columns.")
    return p


def main():
    sys.exit(collect(build_parser().parse_args()))


if __name__ == "__main__":
    main()
