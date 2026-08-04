# Written by Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
Unit test for the Zeng2026_DELSSOME pipelines.

Each test runs one of the bundled example scripts under ``../examples/`` end-to-end
and asserts a zero exit code. The example scripts themselves invoke their
``*check*_results.py`` checkers, which raise if the freshly produced outputs deviate
from the bundled ``reference_output/`` beyond tolerance. So this file only decides
WHICH examples run and WHERE their heavy data comes from; the numeric assertions
live with the examples.

Data resolution: the pFIC and Hopf group-level examples and the GAMLSS example run
from small inputs committed under ``../examples/`` and need nothing else. The
predictor-training example, the group-level DELSSOME example and the individual-level
E/I example additionally need the larger training matrices / pretrained predictors,
which may or may not be committed depending on which copy of the repository this is.
Those are resolved on demand:

  1. if ``../data`` / ``../pretrained_models`` already exist, they are used as-is
     (the case in the stand-alone public repository, where they are committed);
  2. otherwise, if ``$CBIG_TESTDATA_DIR`` is set, the corresponding subtree under
     ``$CBIG_TESTDATA_DIR/stable_projects/fMRI_dynamics/Zeng2026_DELSSOME/`` is
     symlinked into place for the duration of the test and removed afterwards;
  3. otherwise the test is skipped with an explanatory message.

Prerequisites: the ``CBIG_DELSSOME`` conda environment for the group-level tests and
``CBIG_DELSSOME_indiv`` for the individual-level ones. On the CBIG cluster this is
run during release testing, where ``$CBIG_TESTDATA_DIR`` is populated.
"""

import os
import subprocess
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
GROUP_EXAMPLE_SCRIPTS = os.path.join(REPO, "examples", "group_level", "scripts")
INDIV_EXAMPLES = os.path.join(REPO, "examples", "indiv_level")
PROJ_REL = os.path.join("stable_projects", "fMRI_dynamics", "Zeng2026_DELSSOME")
TESTDATA = os.environ.get("CBIG_TESTDATA_DIR")


class TestDELSSOMEExamples(unittest.TestCase):

    def setUp(self):
        self._staged = []

    def tearDown(self):
        for link in self._staged:
            try:
                if os.path.islink(link):
                    os.unlink(link)
            except OSError:
                pass

    def _ensure_staged(self, subdir):
        """Make ``REPO/<subdir>`` available; skip the test if it cannot be located.

        In-repo data is used as-is; otherwise it is symlinked from
        ``$CBIG_TESTDATA_DIR/<PROJ_REL>/<subdir>`` and cleaned up in tearDown.
        Only symlinks are removed, so real in-repo directories are never clobbered.
        """
        target = os.path.join(REPO, subdir)
        if os.path.exists(target):
            return
        if not TESTDATA:
            self.skipTest(f"'{subdir}' is not in the repo and $CBIG_TESTDATA_DIR is unset. "
                          f"Set it and place the data under $CBIG_TESTDATA_DIR/{PROJ_REL}/{subdir}, "
                          "or obtain the data from the stand-alone public repository.")
        src = os.path.join(TESTDATA, PROJ_REL, subdir)
        if not (os.path.isdir(src) or os.path.isfile(src)):
            self.skipTest(f"Test data not found: {src}")
        os.symlink(src, target)
        self._staged.append(target)

    def _run_example(self, path):
        self.assertTrue(os.path.exists(path), f"Missing example script: {path}")
        result = subprocess.run(["bash", path], capture_output=True, text=True)
        self.assertEqual(result.returncode,
                         0,
                         msg=f"{os.path.basename(path)} exited with {result.returncode}\n"
                         f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}")

    def _run_group_example(self, script):
        self._run_example(os.path.join(GROUP_EXAMPLE_SCRIPTS, script))

    def _run_indiv_example(self, script):
        self._run_example(os.path.join(INDIV_EXAMPLES, script))

    # ----- group-level ---------------------------------------------------- #
    def test_cmaes_pfic(self):
        self._run_group_example("CBIG_DELSSOME_CMAES_pFIC_example.sh")

    def test_cmaes_hopf(self):
        self._run_group_example("CBIG_DELSSOME_CMAES_Hopf_example.sh")

    def test_cmaes_pfic_dl(self):
        self._ensure_staged("pretrained_models")
        self._run_group_example("CBIG_DELSSOME_CMAES_pFIC_DL_example.sh")

    def test_predictor_train(self):
        self._ensure_staged("data")
        self._run_group_example("CBIG_DELSSOME_predictor_train_example.sh")

    # ----- individual-level ----------------------------------------------- #
    def test_indiv_estimate_ei(self):
        # The Euler stages need no predictor, but the example also does a
        # DELSSOME-accelerated run when pretrained_models/ holds the
        # individual-level checkpoint -- and reference_output/delssome/ covers
        # that run, so the checker fails if it was skipped. Stage the checkpoint.
        self._ensure_staged("pretrained_models")
        self._run_indiv_example("CBIG_DELSSOME_indiv_EI_example.sh")

    def test_indiv_gamlss(self):
        # Exits 0 with a [skip] notice when Rscript is not on PATH.
        self._run_indiv_example("CBIG_DELSSOME_indiv_gamlss_example.sh")


if __name__ == "__main__":
    unittest.main()
