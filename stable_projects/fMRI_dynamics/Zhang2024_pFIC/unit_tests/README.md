## parameteric Feedback Inhibition Control (pFIC) unit test

The unit tests is for  **CBIG lab only**.

----

## Set up
Please setup the conda environment before running the unit test.
```
cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests
conda env create -f CBIG_pFIC_python_env.yml
conda activate CBIG_pFIC
```

----

## Run

```
python ./unit_tests/test_CBIG_pFIC_unit_test.py
```

The test submits a GPU job (via `CBIG_pbsubmit`, to the `gpuQ-L40S` queue) that runs training + matrix transpose, waits for it to finish, then compares `unit_tests/output/training/training_1.csv` against `unit_tests/ref_results/training/training_1.csv` (must match within `1e-6` summed absolute difference).

**[IMPORTANT]** The reference outputs were generated using an RTX3090 GPU. Jobs now run on `gpuQ-L40S` (the old `gpuserver1/2/3` RTX3090 nodes and `gpuQ` queue no longer exist on the cluster), so results may differ from the reference outputs by more than the comparison tolerance purely due to the different GPU architecture. If the test fails only on this numeric comparison, the reference outputs likely need to be regenerated on the current hardware. You can check which GPU a job used from `unit_tests/output/training/training_iteration_1.txt`.
