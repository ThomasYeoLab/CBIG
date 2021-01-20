# Unit Tests
The unit tests will perform the following steps:
1. Step 1: Converting RSFC data to "documents"
	* Step 1A: Converting RSFC data to "documents" for latent factor estimation
	* Step 1B: Converting RSFC data to "documents" for inferring factor compositions of new participants
2. Step 2: polarLDA
	* Step 2A: Estimating model parameters (or loosely speaking, estimating latent factors)
	* Step 2B: Visualizing latent factors
	* Step 2C: Inferring factor compositions of new participants

NOTE: All the input data and directories in the unit tests **only work for CBIG lab**.

----
## Reference
Siyi Tang*, Nanbo Sun*, Dorothea L. Floris, Xiuming Zhang, Adriana Di Martino, B.T. Thomas Yeo. [Reconciling Dimensional and Categorical Models of Autism Heterogeneity: a Brain Connectomics & Behavioral Study](https://doi.org/10.1016/j.biopsych.2019.11.009). Biological psychiatry, in press.

----
## Data
Data for unit tests are located at `$CBIG_TESTDATA_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/data`

----
## Run Unit Tests & Check Results
NOTE: In the unit tests code, we directly use the copy of executive of polarLDA in `../step2_polarLDA` which was compiled on our CBIG lab server. However, if you are using a different system, we suggest you re-compile the source code of polarLDA by the following command on terminal:
```bash
cp -aR ${CBIG_CODE_DIR}/external_projects/polarlda-c-dist <your_code_dir> # copy to your own code directory
cd <your_code_dir>/polarlda-c-dist/code
make # compile
```
To run the unit tests, go to the unit tests directory `$CBIG_CODE_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/unit_tests/scripts` and run `runtests('CBIG_ASDf_unit_test.m')` in MATLAB.

This MATLAB test script will call the bash script `CBIG_ASDf_unit_test.sh`, which will submit a job to CBIG cluster, and run the steps mentioned above. For Step 2A factor estimation, it will submit another 5 jobs to our cluster to run 5 random initializations for two-factor estimate. After the 5 jobs have finished, the bash script will proceed to Step 2B. Lastly, the test script will compare the outputs with the reference results located at `$CBIG_TESTDATA_DIR/stable_projects/disorder_subtypes/Tang2019_ASDFactors/results`.

----
## Bugs and Questions
Please contact Siyi Tang at tangsy935@gmail.com, Nanbo Sun at sun464879934@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
