# LDA -- General

**Note** These are general functions that can be applied to new data. To replicate our PNAS results, learn the overall procedure, or see examples of how to call these functions, please refer to folder `../replicatePNAS/`.

---
## LDA Source Code

A copy of the LDA source code, written by David Blei, is in `${CBIG_CODE_DIR}/external_packages/lda-c-dist/`. To keep that copy intact, the wrapper `../replicatePNAS/CBIG_LDA_wrapper.m` will copy that folder and compile it here. Therefore, when you run the wrapper, there will be a folder named `lda-c-dist` generated here, containing both the source code and executive.

----
## Step 1: Converting VBM Results into Documents

Relevant File: `CBIG_brain2doc.m`

This is the first step that converts the VBM results (GM density images) into "brain documents", which LDA can analyze.

Specifically, it will first take log10 of the voxel values, regress out nuisance variables for all images w.r.t. the "reference group", z-normalize each voxel of all images w.r.t. the same reference group, and finally output the document files, which will later serve as the input of LDA.

To see detailed usage, type `help CBIG_brain2doc` in MATLAB. For an example of calling it, see `../replicatePNAS/CBIG_LDA_wrapper.m`.

----
## Step 2: Estimating Model Parameters (Atrophy Factors)

Relevant Files: `CBIG_LDA_est.sh`, `CBIG_LDA_est_job.sh`, `CBIG_LDA_infSettings.txt`, `CBIG_waitUntilFinished.m`

In the second step, we estimate the LDA model parameters (one of which is Pr(Voxel | Factor), probabilistic maps of atrophy factors) from the selected group of images (now documents).

`CBIG_LDA_est.sh` is the master script that either submits `CBIG_LDA_est_job.sh` to your cluster or runs it serially. `CBIG_LDA_est_job.sh` runs the LDA executive from `./lda-c-dist/` and uses `CBIG_LDA_infSettings.txt` as inference parameters. Finally, `CBIG_waitUntilFinished.m` is just a helper function that holds off the next steps until all the random initializations are done.

Since LDA cannot guarantee global optimum, we need to run multiple (e.g., 20) random initializations and finally pick the best one. For more technical details on LDA, please refer to the original paper
>Blei DM, Ng AY, Jordan MI (2003) Latent Dirichlet allocation. *J Machine Learning Res* 3:993-1022.

For an example of calling `CBIG_LDA_est.sh`, see `../replicatePNAS/CBIG_LDA_wrapper.m`.

----
## Step 3: Visualizing Estimated Atrophy Factors

Relevant File: `CBIG_visualizeFactors.m`

This function selects the best random initilization (that gives the highest likelihood) and visualizes the estimated factors as probabilistic atrophy maps. This visualization involves writing probability values into the brain. the In addition, it also normalizes the gamma in LDA to produce factor composition that sum to 1 (e.g., `[0.7, 0.2, 0.1]`).

To see detailed usage, type `help CBIG_visualizeFactors` in MATLAB. For an example of calling it, see `../replicatePNAS/CBIG_LDA_wrapper.m`.

----
## Step 4: Inferring New subjects with Learned Model

Relevant File: `CBIG_LDA_inf.sh`

This function infers factor compositions of new subjects with the learned LDA model. In other words, this is the function to use when you have estimated your factors and want to extract the factor composition of a new subject with you model. For example in our PNAS paper, we estimated our LDA model parameters (or loosely speaking, the atrophy factors) with 188 ADNI-1 AD patients and then used this script to infer factor compositions of 190 amyloid-positive nondemented participants (unseen subjects).

For an example of calling it, see `../replicatePNAS/CBIG_LDA_wrapper.m`.

----
## How do lines in the factor composition output file correspond to the images?

The principle here is that lines in the factor composition output file always match with order of the input documents (images). Hence, if you input a filename list to `../lib/CBIG_brain2doc.m`, lines in the factor composition output file will correspond to that list. That is, the first line of factor compositions will correspond to the first image listed in the filename list. Otherwise, they will match the order of the `concatOrder` file.
