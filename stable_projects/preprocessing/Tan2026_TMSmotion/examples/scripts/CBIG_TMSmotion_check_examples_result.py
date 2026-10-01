'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

This script is called by examples/CBIG_TMSmotion_examples_wrapper.sh (not meant
to be run directly by the user), e.g.:

Examples:

bash $TAN2026_TMSMOTION_DIR/examples/CBIG_TMSmotion_examples_wrapper.sh
'''

import os
import numpy as np


def check_npy(output_file, ref_file, rtol=1e-5, atol=1e-8):
    '''Compare two saved numpy arrays for numerical closeness.

    Loads `output_file` and `ref_file` as numpy arrays, verifies both
    exist, checks that their shapes match, and asserts their values are
    close within the given tolerances. Prints a confirmation message on
    success.

    Args:
        output_file (str): Path to the `.npy` file containing the
            generated output.
        ref_file (str): Path to the `.npy` file containing the reference
            output.
        rtol (float, optional): Relative tolerance passed to
            `np.allclose`. Defaults to 1e-5.
        atol (float, optional): Absolute tolerance passed to
            `np.allclose`. Defaults to 1e-8.

    Raises:
        FileNotFoundError: If `output_file` or `ref_file` does not
            exist.
        ValueError: If the arrays' shapes differ, or their values are
            not close within tolerance.
    '''

    if not os.path.exists(output_file):
        raise FileNotFoundError(f"Missing output file: {output_file}")

    if not os.path.exists(ref_file):
        raise FileNotFoundError(f"Missing reference file: {ref_file}")

    arr_out = np.load(output_file)
    arr_ref = np.load(ref_file)

    if arr_out.shape != arr_ref.shape:
        raise ValueError(f"Shape mismatch:\n"
                         f"Output shape: {arr_out.shape}\n"
                         f"Ref shape:    {arr_ref.shape}")

    if not np.allclose(arr_out, arr_ref, rtol=rtol, atol=atol):
        diff = np.max(np.abs(arr_out - arr_ref))
        raise ValueError(f"Arrays differ! Max absolute difference = {diff}")

    print(f"✓ PASS: {os.path.basename(output_file)}")


def main():
    '''Validate example outputs against reference outputs for all subjects.

    For each subject in `["sub-01", "sub-02"]`, compares the generated
    `.npy` output in `<TAN2026_TMSMOTION_DIR>/examples/output/` against
    the reference `.npy` file in
    `<TAN2026_TMSMOTION_DIR>/examples/ref_output/` via `check_npy`.
    Prints a success message if all comparisons pass.

    Raises:
        EnvironmentError: If the `TAN2026_TMSMOTION_DIR` environment
            variable is not set.
    '''

    TAN2026_TMSMOTION_DIR = os.environ.get("TAN2026_TMSMOTION_DIR")
    if TAN2026_TMSMOTION_DIR is None:
        raise EnvironmentError("Please set TAN2026_TMSMOTION_DIR.")

    output_dir = os.path.join(TAN2026_TMSMOTION_DIR, "examples", "output")

    ref_dir = os.path.join(TAN2026_TMSMOTION_DIR, "examples", "ref_output")

    subjects = ["sub-01", "sub-02"]

    for subID in subjects:
        output_file = os.path.join(output_dir, f"{subID}.npy")
        ref_file = os.path.join(ref_dir, f"{subID}.npy")

        check_npy(output_file, ref_file)

    print("All example outputs match reference outputs.")
    print("Examples validation completed successfully.\n")


if __name__ == "__main__":
    main()
