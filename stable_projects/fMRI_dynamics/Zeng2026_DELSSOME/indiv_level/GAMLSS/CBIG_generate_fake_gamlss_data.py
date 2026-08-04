#!/usr/bin/env python3
# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Synthesise a small tabular dataset that has the same column layout the
``CBIG_fit_gamlss_model.R`` script expects, so the GAMLSS pipeline can be
run without needing access to the real 12,005-participant E/I-ratio table.

Required columns (consumed by CBIG_gamlss_utils.R):

* ``age_in_years``  — numeric age.
* ``numeric_sex``    — 0 (male) or 1 (female) — used as a fixed effect.
* ``site_scanner``   — categorical site x scanner — modelled as a random effect.
* ``ei_ratio_mean``  — the response variable being modelled (the example uses
                       a non-linear lifespan trajectory + sex offset + noise).

Optionally:
* ``subject_id``, ``sex``, ``site``, ``scanner``  — kept for traceability.
* ``ei_ratio_roi_{0..n_roi-1}``  — fake per-ROI E/I ratios, in case the user
                                   wants to also fit GAMLSS to each ROI.

Output is written to ``--out`` (default ``fake_gamlss_input.csv`` next to the
R scripts).
"""

import argparse
import os
import numpy as np
import pandas as pd


def lifespan_ei_curve(age, sex):
    """Mock normative E/I curve: decrease through development + late-life rise."""
    base = (
        2.6
        - 0.025 * np.maximum(age - 5.0, 0.0)
        + 0.0006 * np.maximum(age - 50.0, 0.0) ** 1.5
    )
    base = base + 0.05 * sex  # higher in females (as in the paper)
    return base


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        default=os.path.join(os.path.dirname(__file__), "fake_gamlss_input.csv"),
        help="Output CSV path.",
    )
    parser.add_argument("--n", type=int, default=600,
                        help="Number of fake participants.")
    parser.add_argument("--n-roi", type=int, default=68)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    N = args.n

    # Age distribution that roughly mimics the paper (4.5 to 98.4 yrs).
    age = rng.uniform(4.5, 98.4, size=N)
    sex_label = rng.choice(["male", "female"], size=N)
    numeric_sex = (sex_label == "female").astype(int)

    # Half a dozen sites x scanners to exercise the random-effect part.
    sites = rng.choice(
        ["WashU", "MGH", "UPenn", "NUS", "NTU", "LMU"], size=N
    )
    scanners = rng.choice(["Siemens", "GE"], size=N)
    site_scanner = np.array([f"{s}_{sc}" for s, sc in zip(sites, scanners)])

    # Add a small per-site offset so harmonisation has a job to do.
    site_offsets = {s: rng.normal(0.0, 0.05) for s in np.unique(site_scanner)}
    site_eff = np.array([site_offsets[ss] for ss in site_scanner])

    noise = rng.normal(0.0, 0.05, size=N)
    ei_ratio_mean = lifespan_ei_curve(age, numeric_sex) + site_eff + noise

    df = pd.DataFrame({
        "subject_id": [f"fake-{i:04d}" for i in range(N)],
        "age_in_years": age,
        "sex": sex_label,
        "site": sites,
        "scanner": scanners,
        "site_scanner": site_scanner,
        "numeric_sex": numeric_sex,
        "ei_ratio_mean": ei_ratio_mean,
    })

    # Per-ROI fake EI ratios = mean +/- per-ROI offsets.
    roi_offsets = rng.normal(0.0, 0.1, size=args.n_roi)
    for r in range(args.n_roi):
        df[f"ei_ratio_roi_{r}"] = ei_ratio_mean + roi_offsets[r] + rng.normal(
            0.0, 0.03, size=N
        )

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"Wrote {N} fake participants to {args.out}")


if __name__ == "__main__":
    main()
