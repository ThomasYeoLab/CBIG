# Parametric Mean Field Model Analysis Main Results
* `input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `scripts` this folder contains the scripts to generate the main results of the paper. 
    * `CBIG_pMFM_basic_functions_main.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.
    * `CBIG_pMFM_step1_training_main.py` this file is the wrapper to perform parameter estimation process. The process performs 10 random initializations and each random initialization takes 500 iterations. In total, there are 5000 candidate parameter sets. The estimated time for one random initialization is 2 days. Total time used for the whole training process (step1) is around 20 days.
    * `CBIG_pMFM_step2_validation_main.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_main.py`. The 5000 candidate parameter sets are evaluated in the validation set to obtain top 10 candidate parameter sets.
    * `CBIG_pMFM_step3_test_main.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_main.py`. The script generates 1000 simulations for each parameter set and computes the averaged cost.
    * `CBIG_pMFM_step4_generate_simulated_fc_fcd.py` this file is the wrapper to generate simulated FC and FCD based on MFM and estimated model parameters from `CBIG_pMFM_step3_test_main.py`
    * `CBIG_pMFM_step5_generate_STDFCD_correlation_main.m` this file is the function to compute the empirical and simulated SWSTD-FCD correlation
    * `CBIG_pMFM_step6_fitmodel.py` this file is the function to compute the distribution of FCD mean and the threshold if FCD mean follows the bi-modal distribution
    * `CBIG_pMFM_step6_SWSTD_state_main.m` this file is the function to perform SWSTD state analysis which corresponds to figure 5D and 5E in the paper.
    * `CBIG_pMFM_step7_perturbation_analysis.py` this file is the wrapper to perform perturbation analysis. This step takes about 1.5 days to finish.
    * `CBIG_pMFM_step8_gene_expression_analysis_desikan.m` this file is the wrapper to perform gene expression analysis 

# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part1_pMFM_main/scripts` and follow the step order from step 1 to step 7 to run the scripts.
* For Python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts. 
* For Matlab scripts, directly run the scripts with proper inputs.

# Results
* Top 10 model parameters and their corresponding costs are store in a CSV file called `test_all.csv` in the `output/step3_test_results` folder.
    * There are 10 columns in `test_all.csv` file. The first column contains the parameters and costs for the top 1 result; second column contains the second top result and so on.
    * The first row represents the random initializaiton index.
    * The second row represents the index of generating in the training process
    * 3-5 rows represents the training FC costs, training FCD KS statistics and training total costs
    * 6-8 rows represents the validation FC costs, validation FCD KS statistics and validation total costs
    * 9-11 rows represents the test FC costs, test FCD KS statistics and test total costs
    * 12-79 rows represents recurrent connection w of 68 ROIs
    * 80-147 rows represents external input I of 68 ROIs
    * 148 row represents global amplitude G
    * 149-216 rows represents noise amplitude sigma of 68 ROIs
