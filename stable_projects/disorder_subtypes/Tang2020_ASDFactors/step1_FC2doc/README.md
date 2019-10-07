# Converting RSFC data to "documents"
This step converts the whole-brain RSFC data into "documents", which can be analyzed by the Bayesian model proposed in our paper (we call it "polarLDA").

Specifically, this step regresses out nuisance variables from all participants' RSFC data and z-normalizes with respect to the control participants, and finally writes into document files that will later serve as the input of polarLDA model.

----
## "Documents" for Estimating Latent Factors
Relevant files:`CBIG_ASDf_FC2doc.m`,`CBIG_ASDf_FC2doc_estFactors_wrapper.m`

To prepare documents for latent factor estimation, the function `CBIG_ASDf_FC2doc.m` performs regression of nuisance variables (e.g., age, sex, motion & sites) with a GLM estimated from only the control participants. It then performs z-normalization with respect to controls, times 10 of the z-scores and discretizes the values.

`CBIG_ASDf_FC2doc_estFactors_wrapper.m` is the wrapper function to call `CBIG_ASDf_FC2doc.m` to obtain the z-normalized, discretized RSFC data of ABIDE-II+GENDAAR participants.

----
## "Documents" for Inferring Factor Compositions of New Participants
Relevant files:`CBIG_ASDf_FC2doc_forInference.m`,`CBIG_ASDf_FC2doc_infFactorComp_wrapper.m`

To prepare documents for inferring factor compositions of new participants, the function `CBIG_ASDf_FC2doc_forInference.m` performs nuisance variable regression with respect to control participants, and then z-normalization with respect to controls in the "reference sample", i.e., the sample used to estimate the latent factors. It also times 10 of the z-scores and discretizes the values.

Note that we perform nuisance variable regression with respect to control participants in one sample (the sample for factor composition inference), but z-normalization with respect to controls in the "reference sample". This is because nuisance variables include acquisition sites, and different samples might have different acquisition sites. Therefore, regression is not performed with respect to the "reference sample".

`CBIG_ASDf_FC2doc_infFactorComp_wrapper.m` is the wrapper function to call `CBIG_ASDf_FC2doc_forInference.m` to obtain the z-normalized, discretized RSFC data of ABIDE-I participants (z-normalization w.r.t. ABIDE-II+GENDAAR control participants).

