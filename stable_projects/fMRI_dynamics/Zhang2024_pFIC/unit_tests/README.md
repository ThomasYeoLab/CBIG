## parameteric Feedback Inhibition Control (pFIC) unit test

The unit tests is for  **CBIG lab only**.

Run the following commands to get started

```
ssh headnode

cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC

bash unit_tests/CBIG_pFIC_unit_test_submission.sh
```
or in Matlab, `cd` to the `unit_tests` folder, then run

```
runtests(pwd)
```
**[IMPORTANT]** Note that the reference outputs are generated using RTX3090, make sure that the job is submitted to the node  `gpuserver1`, `gpuserver2` or `gpuserver3`. Submission to the node `gpuserver4` will lead to failing the unit test as it is equipped with a different GPU card. When the job starts running, you can check which GPU it's using from `examples/output/training/training_iteration_1.txt` .
