# Written by Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# Co-authored by Fang Tian

"""Regression check for the individual-level E/I example.

Diffs every CSV under ``reference_output/euler/<sub_id>/`` against the matching
file under ``output/euler/<sub_id>/``, and does the same for ``delssome/`` when
a reference for it exists. The checked cohort is defined by whatever the
reference snapshot contains, so adding a subject to the reference automatically
adds it to the check. ``reference_output/gamlss/`` is deliberately out of scope
-- CBIG_DELSSOME_indiv_check_gamlss_example_results.py owns it.

Tolerance is 1e-6 on the sum of absolute differences, i.e. effectively
bit-exact. CMA-ES is seeded, so a mismatch means either a code change or a
different BLAS / CPU dispatch path -- regenerate the references on the reference
machine rather than loosening this.
"""

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_ROOT = os.path.join(HERE, "output")
REFERENCE_ROOT = os.path.join(HERE, "reference_output")
# The subtrees this checker owns, one per simulation engine the demo runs.
ENGINE_SUBDIRS = ("euler", "delssome")
TOL = 1e-6


def compare_csv(ref_path, out_path):
    """Return None when the two CSVs agree, else a one-line reason."""
    if not os.path.isfile(out_path):
        return f"missing output {out_path}"
    ref = pd.read_csv(ref_path)
    out = pd.read_csv(out_path)
    if list(ref.columns) != list(out.columns):
        return (f"columns differ: reference {list(ref.columns)} vs "
                f"output {list(out.columns)}")
    if len(ref) != len(out):
        return f"row count differs: reference {len(ref)} vs output {len(out)}"
    numeric = ref.select_dtypes(include=[np.number]).columns
    diff = float(np.abs(out[numeric].to_numpy(dtype=float)
                        - ref[numeric].to_numpy(dtype=float)).sum())
    if not np.isfinite(diff):
        return "output contains NaN or inf"
    if diff > TOL:
        return f"sum |difference| = {diff:.3e} > {TOL:.0e}"
    return None


def collect_reference_csvs():
    """``[(relative_path, ref_path), ...]`` for every reference CSV this demo owns.

    Only the ``euler/`` and ``delssome/`` subtrees are walked. The rest of
    ``reference_output/`` belongs to other demos -- notably ``gamlss/``, which is
    produced and checked by CBIG_DELSSOME_indiv_gamlss_example.sh -- and claiming
    it here would fail whenever that demo has not run.
    """
    found = []
    for subdir in ENGINE_SUBDIRS:
        for root, _dirs, files in os.walk(os.path.join(REFERENCE_ROOT, subdir)):
            for name in sorted(files):
                if name.endswith(".csv"):
                    ref_path = os.path.join(root, name)
                    found.append((os.path.relpath(ref_path, REFERENCE_ROOT), ref_path))
    return sorted(found)


NO_REFERENCE_HINT = """\
[check_EI] No reference CSVs under
    {root}

The demo itself succeeded and its results are under output/euler/. Only the
comparison against the released reference snapshot could not run, because that
snapshot is missing from this copy of the repository."""


def main():
    references = collect_reference_csvs() if os.path.isdir(REFERENCE_ROOT) else []
    if not references:
        sys.exit(NO_REFERENCE_HINT.format(root=REFERENCE_ROOT))

    failures = {}
    for rel_path, ref_path in references:
        reason = compare_csv(ref_path, os.path.join(OUTPUT_ROOT, rel_path))
        if reason is not None:
            failures[rel_path] = reason

    if failures:
        for rel_path, reason in failures.items():
            print(f"  {rel_path}: {reason}")
        raise Exception(
            f"Failed. E/I example results differ from the reference for "
            f"{len(failures)}/{len(references)} CSV(s).")
    print(f"Passed! E/I example results match the reference for all "
          f"{len(references)} CSV(s).")


if __name__ == "__main__":
    main()
