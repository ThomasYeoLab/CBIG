# Examples

## Background
These examples run code used to simulate rs-fMRI runs for 2 test participants.

## Data
Input data is stored under `input` folder.
* `window_list_1.txt` contains data on each 1-minute non-overlapping rs-fMRI window.
  - Each row in this file corresponds to 1 window from 1 participant.
  - This file contains data for 2 example participants.
* `pipeline.json` contains example simulation parameters (i.e., # sampling iterations).
  - For examples, # sampling iterations were reduced to speed up computation.

## Run
* Please follow all steps under "Setup" section from `${TAN2026_TMSMOTION_DIR}/README.md`.
* Run on terminal:
    ```
    bash ${TAN2026_TMSMOTION_DIR}/examples/CBIG_TMSmotion_examples_wrapper.sh
    ```

Outputs will be generated under `$TAN2026_TMSMOTION_DIR/examples/output`.
