This folder contains wrapper script `scripts/CBIG_ASDf_replication.sh` to replicate main results in our paper.

The script will perform the following steps:
1. Step 1: Converting RSFC data to "documents"
	* Step 1A: Converting RSFC data to "documents" for latent factor estimation
	* Step 1B: Converting RSFC data to "documents" for inferring factor compositions of new participants
2. Step 2: polarLDA
	* Step 2A: Estimating model parameters (or loosely speaking, estimating latent factors)
	* Step 2B: Visualizing latent factors
	* Step 2C: Inferring factor compositions of new participants

NOTE: All the input data and directories in this folder **only work for CBIG lab**. Due to MATLAB version difference, the results generated using this script might not be exactly the same as the initial results in our paper, but the difference should be minor.

----
## Reference
Siyi Tang*, Nanbo Sun*, Dorothea L. Floris, Xiuming Zhang, Adriana Di Martino, B.T. Thomas Yeo. [Reconciling Dimensional and Categorical Models of Autism Heterogeneity: a Brain Connectomics & Behavioral Study](https://doi.org/10.1016/j.biopsych.2019.11.009). Biological psychiatry, in press.

----
## Data
Data used in this folder are located at `$CBIG_REPDATA_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/data`
### Subjects list
You can find the subjects list in the following files:
1. ABIDE-II+GENDAAR: `$CBIG_REPDATA_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/data/subInfo_654.csv`
2. ABIDE-I: `$CBIG_REPDATA_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/data/subInfo_316.csv`

----
## Run Replication Script
On terminal, specify the output directory and call the script `$CBIG_CODE_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/replication/scripts/CBIG_ASDf_replication.sh <your_output_dir>`.

This script will submit a job to `circ-spool` cluster, and run the steps mentioned above. For Step 2A factor estimation, it will submit another 100 jobs to our cluster to run 100 random initializations for two-factor estimate, as well as another 100 jobs to run 100 random initializations for three-factor estimate. After the 200 factor estimation jobs have finished, the script will proceed to Step 2B. 

After the jobs have finished, you may compare your outputs with the reference results located at `$CBIG_REPDATA_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors/results`.

----
## Bugs and Questions
Please contact Siyi Tang at tangsy935@gmail.com, Nanbo Sun at sun464879934@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
