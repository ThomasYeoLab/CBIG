# Estimating Latent ASD Factors, Visualizing Factors & Inferring Factor Compositions of New Participants
This step applies the Bayesian model proposed in our paper (we call it "polarLDA model") to estimate latent ASD factors. Moreover, this step also allows you to visualize the factor-specific RSFC patterns of the latent factors (i.e., E(RSFC patterns|Factor)). Lastly, you can infer factor compositions of new participants with the estimated model parameters.

For examples of using the functions in this folder, please see `../examples/scripts/CBIG_ASDf_example_script.sh` and `../unit_tests/scripts/CBIG_ASDf_unit_test.sh`.

----

## Source Code of polarLDA Model
The source code of our polarLDA model has been released in `$CBIG_CODE_DIR/external_packages/polarlda-c-dist/code`. For more technical details on polarLDA, please refer to supplemental materials in our paper and `$CBIG_CODE_DIR/external_packages/polarlda-c-dist/doc/polarLDA.pdf`.

----

## Step A: Estimating Model Parameters (Latent Factors)
Relevant files: `CBIG_ASDf_polarLDA_est.csh`, `CBIG_ASDf_polarLDA_est_job.sh`, `CBIG_ASDf_polarLDA_infSettings.txt`,`polarLDA`

In step A, we estimate the model parameters (or loosely speaking, the latent factors) from ASD participants' z-normalized, discretized RSFC data obtained from `../step1_FC2doc`. 

The main script is `CBIG_ASDf_polarLDA_est.csh`. It runs the polarLDA executive `polarLDA` using `CBIG_ASDf_infSettings.txt` as inference settings. Because the algorithm does not guarantee global optimal, you need to run multiple random initializations (e.g., 100) and pick a solution as the final estimate. If you have a cluster, `CBIG_ASDf_polarLDA_est.csh` will call the script `CBIG_ASDf_polarLDA_est_job.sh` and submit jobs to your cluster, so that you can run different random initialization in parallel. Otherwise, different random initializations will run in series.

The file `polarLDA` is the executive of polarLDA compiled on our CBIG lab server. If you are using a different system, we suggest you re-compile the source code of polarLDA by the following commands on terminal:
```bash
# copy to your own code directory to keep the common space clean
cp -aR ${CBIG_CODE_DIR}/external_projects/polarlda-c-dist <your_code_dir>
# compile
cd <your_code_dir>/polarlda-c-dist/code; make
```

----

## Step B: Visualizing Latent Factors
Relevant files: `CBIG_ASDf_computeCorrBetweenRuns.m`, `CBIG_ASDf_corrTwoRuns.m`,`CBIG_ASDf_gamma2table.m`,`CBIG_ASDf_mean2mat.m`,`CBIG_ASDf_hunMatch.m`,`CBIG_ASDf_plotCorrWithBest.m`,`CBIG_ASDf_plotLogLike.m`,`CBIG_ASDf_visualizeFactors.m`

In step B, the main function is `CBIG_ASDf_visualizeFactors.m`. This function first computes the pairwise correlation between all the solutions (from different random initializations), and selects the solution with the highest average correlation with all other solutions as the final estimate. Next, it visualizes the estimated latent factors by computing the expected factor-specific hypo/hyper RSFC patterns (i.e., E(RSFC patterns|Factor)), and plotting each factor as a 419 x 419 matrix. In addition, it also normalizes the estimated gamma from polarLDA model so that each participant's factor composition sums up to 1 (e.g., [0.7, 0.2, 0.1]).

----

## Step C: Inferring Factor Compositions of New Participants
Relevant files: `CBIG_ASDf_gamma2table.m`,`CBIG_ASDf_polarLDA_inf.csh`,`CBIG_ASDf_polarLDA_inference_wrapper.sh`

In step C, script `CBIG_ASDf_polarLDA_inf.csh` infers factor compositions of new participants with the learned polarLDA model. In other words, this is the function to use if you have estimated your latent factors and want to extract the factor compositions of a new participant. For example, in our paper, we estimated our polarLDA model parameters with 306 ASD participants in ABIDE-II+GENDAAR datasets, and inferred factor compositions of 166 ASD participants (new participants) in ABIDE-I datasets.

The wrapper script `CBIG_ASDf_polarLDA_inference_wrapper.sh` calls `CBIG_ASDf_polarLDA_inf.csh` to infer gamma (i.e., unnormalized factor compositions) using the final estimated model parameters in previous steps. It also normalizes gamma so that the factor composition of each participant sum up to 1.
