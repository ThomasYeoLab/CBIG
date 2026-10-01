'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

This script runs unit tests to check the simulations of rs-fMRI runs
for two test participants.

Examples:
    cd $TAN2026_TMSMOTION_DIR/unit_tests
    source CBIG_init_conda; conda activate CBIG_TMSmotion
    python test_CBIG_TMSmotion_unit_test.py
'''

import os
import unittest
import subprocess
import numpy as np

# ---------------------------------------------------------------------
# Environment setup
# ---------------------------------------------------------------------
UNIT_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(UNIT_TESTS_DIR)
EXAMPLES_SCRIPT = os.path.join(PROJECT_DIR, "examples", "scripts",
                               "CBIG_TMSmotion_examples.py")
OUTPUT_DIR = os.path.join(UNIT_TESTS_DIR, "output")
REF_DIR = os.path.join(PROJECT_DIR, "examples", "ref_output")

SUBJECTS = ["sub-01", "sub-02"]


# ---------------------------------------------------------------------
# Test class
# ---------------------------------------------------------------------
class TestTMSmotion(unittest.TestCase):

    def test_simulation(self):
        """Verify that the examples simulation script reproduces reference output.

        Removes any pre-existing per-subject `.npy` output files in
        `OUTPUT_DIR`, runs `CBIG_TMSmotion_examples.py` (invoked with the
        "unit_tests" output directory name) as a subprocess, and asserts
        it exits successfully. For each subject in `SUBJECTS`, asserts
        that both the newly generated output file and the corresponding
        reference file in `REF_DIR` exist, that their array shapes
        match, and that their values are numerically close (via
        `np.allclose` with `rtol=1e-5`, `atol=1e-8`).

        Raises:
            RuntimeError: If the simulation subprocess exits with a
                non-zero return code.
            AssertionError: If any output file is missing, or its shape
                or values differ from the reference.
        """
        # Clear previous outputs
        for subID in SUBJECTS:
            out_file = os.path.join(OUTPUT_DIR, f"{subID}.npy")
            if os.path.exists(out_file):
                os.remove(out_file)

        # Run the examples simulation script
        cmd = f"python {EXAMPLES_SCRIPT} unit_tests"
        result = subprocess.run(cmd,
                                shell=True,
                                capture_output=True,
                                text=True)
        if result.returncode != 0:
            raise RuntimeError("Simulation script failed:\n"
                               f"STDOUT:\n{result.stdout}\n"
                               f"STDERR:\n{result.stderr}")

        # Compare outputs to reference
        for subID in SUBJECTS:
            out_file = os.path.join(OUTPUT_DIR, f"{subID}.npy")
            ref_file = os.path.join(REF_DIR, f"{subID}.npy")

            self.assertTrue(os.path.exists(out_file),
                            f"Output file missing: {out_file}")
            self.assertTrue(os.path.exists(ref_file),
                            f"Reference file missing: {ref_file}")

            arr_out = np.load(out_file)
            arr_ref = np.load(ref_file)

            self.assertEqual(arr_out.shape,
                             arr_ref.shape,
                             msg=(f"Shape mismatch for {subID}: "
                                  f"{arr_out.shape} vs {arr_ref.shape}"))

            self.assertTrue(
                np.allclose(arr_out, arr_ref, rtol=1e-5, atol=1e-8),
                f"Simulation output differs from reference for {subID}")


if __name__ == "__main__":
    unittest.main()
