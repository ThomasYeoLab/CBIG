# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""Regression check for the GAMLSS example.

Diffs every CSV that ``CBIG_fit_gamlss_model.R`` wrote against the bundled
reference snapshot under ``reference_output/gamlss/<config>/``. Only the
released default config ships a snapshot; with any other ``CONFIG_NAME`` the
comparison is skipped.

Tolerance is 1e-4 absolute. Unlike the E/I check this is deliberately loose:
the GAMLSS RS() iteration is seeded but its convergence path depends on the
BLAS and the R / gamlss build.
"""

import os
import sys

import pandas as pd
from pandas.testing import assert_frame_equal

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG_NAME = "SHASHo2_bs_df4.json"
CONFIG_NAME = os.environ.get("CONFIG_NAME", DEFAULT_CONFIG_NAME)
CONFIG_STEM = CONFIG_NAME.removesuffix(".json")
OUTPUT_ROOT = os.path.join(HERE, "output", "gamlss", CONFIG_STEM)
REFERENCE_ROOT = os.path.join(HERE, "reference_output", "gamlss", CONFIG_STEM)
ATOL = 1e-4


# A non-default CONFIG_NAME is a legitimate thing to run, and no snapshot is
# shipped for one, so that case is a skip rather than a failure.
OTHER_CONFIG_NOTICE = """\
[skip] No reference snapshot for {config}, so the comparison is skipped. Only
       the released default ({default}) ships one. The GAMLSS fit itself
       succeeded; its results are under output/gamlss/{stem}/."""

NO_REFERENCE_HINT = """\
[check_gamlss] No reference CSVs under
    {root}

The GAMLSS fit itself succeeded and its results are under output/gamlss/{stem}/.
Only the comparison against the released reference snapshot could not run,
because that snapshot is missing from this copy of the repository."""


def main():
    ref_csvs = sorted(f for f in os.listdir(REFERENCE_ROOT)
                      if f.endswith(".csv")) if os.path.isdir(REFERENCE_ROOT) else []
    if not ref_csvs:
        if CONFIG_NAME != DEFAULT_CONFIG_NAME:
            print(OTHER_CONFIG_NOTICE.format(config=CONFIG_NAME,
                                             default=DEFAULT_CONFIG_NAME,
                                             stem=CONFIG_STEM))
            return
        sys.exit(NO_REFERENCE_HINT.format(root=REFERENCE_ROOT, stem=CONFIG_STEM))

    failures = {}
    for name in ref_csvs:
        out_path = os.path.join(OUTPUT_ROOT, name)
        if not os.path.isfile(out_path):
            failures[name] = f"missing output {out_path}"
            continue
        try:
            assert_frame_equal(pd.read_csv(out_path),
                               pd.read_csv(os.path.join(REFERENCE_ROOT, name)),
                               atol=ATOL, rtol=0, check_dtype=False)
        except AssertionError as exc:
            failures[name] = str(exc).splitlines()[0]

    if failures:
        for name, reason in failures.items():
            print(f"  {name}: {reason}")
        raise Exception(
            f"Failed. GAMLSS results differ from the reference for "
            f"{len(failures)}/{len(ref_csvs)} CSV(s).")
    print(f"Passed! GAMLSS results match the reference for all "
          f"{len(ref_csvs)} CSV(s).")


if __name__ == "__main__":
    main()
