# Unit tests

## Background
This unit test checks the simulation of rs-fMRI runs for 2 test participants.

## Data
Input data is stored under `$TAN2026_TMSMOTION_DIR/examples/input` folder.
* `window_list_1.txt` contains data on each 1-minute non-overlapping rs-fMRI window.
  - Each row in this file corresponds to 1 window from 1 participant.
  - This file contains data for 2 test participants.
* `pipeline.json` contains example simulation parameters (e.g., # sampling iterations).
  - For unit tests, # sampling iterations were reduced to speed up computation.

## Run
The unit test is for  **CBIG lab only**.

* Please follow all steps under "Setup" section from `${TAN2026_TMSMOTION_DIR}/README.md`.
* Please run on terminal:
```
cd $TAN2026_TMSMOTION_DIR/unit_tests
source CBIG_init_conda; conda activate CBIG_TMSmotion
python test_CBIG_TMSmotion_unit_test.py
```

Outputs will be generated under `$TAN2026_TMSMOTION_DIR/unit_tests/output`.
