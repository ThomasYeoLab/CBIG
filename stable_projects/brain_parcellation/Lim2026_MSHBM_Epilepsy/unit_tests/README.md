# Unit Test: Individual-Specific Epilepsy Network Parcellation

The unit test runs the example wrapper provided in the `examples/` folder and verifies
that the output parcellation exactly matches the reference results.

The unit test will:

+ Run `CBIG_MSHBM_Epilepsy_wrapper()` on the 5 example resting-state fMRI runs and
  motion censor files provided in `examples/input/`. The wrapper calls
  `CBIG_MSHBM_Epilepsy_LI`, which uses `Kong2019_MSHBM` scripts as dependencies.

+ Generate an individual-specific 15-network cortical parcellation using the NIH
  34-subject drug-resistant epilepsy MS-HBM model.

+ Compare the output parcellation labels (`lh_labels`, `rh_labels`) against the
  reference stored in `examples/results/` and assert that the maximum difference is 0.

+ Compare the language laterality index value (`LI_lang`) against the reference stored
  in `examples/results/LI.mat` and assert that the difference is less than 1e-6.

**Dependency:** `Kong2019_MSHBM` must be present under `$CBIG_CODE_DIR/stable_projects/brain_parcellation/`.

Note: the unit test is **for CBIG lab only**.

----

## Bugs and Questions

Please contact Mervyn Lim Jun Rui at mervynlim@u.nus.edu.
