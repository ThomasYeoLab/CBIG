# MMLDA - Generic Functions

This folder contains generic function that can be applied to new data. To replicate our results, learn the overall procedure, or see examples of how to call these function, please refer to folder `../unit_tests/replicate/`. 

----

## MMLDA Source Code
The MMLDA source code is written by C language and modified from LDA code by David Blei. If you want to use it in your machine, we suggest you to recompile it 

```
cd $CBIG_CODE_DIR/external_packages/mmlda-c-dist
make
```
even though we provide the binary executable file (compiled in CBIG lab)

```
MMLDA
```

----

## Step1: Converting VBM results and behavioral scores into documents
1. We need to use reference group to create reference paramters for regression and z-normalization.

```
CBIG_MMLDA_brain_to_zscore_create_refpara.m
CBIG_MMLDA_behavior_to_zscore_create_refpara.m
```

2. Given a new cohort, we can convert VBM results and behavioral scores to zscores by using reference paramters.

```
CBIG_MMLDA_brain_to_zscore.m
CBIG_MMLDA_behavior_to_zscore.m
```

3. Finally, we convert brain and behavior zscore to documents used by MMLDA

```
CBIG_MMLDA_brain_behavior_zscore_to_doc.m
```
To see detailed usage, type `help <function_name>` in MATLAB. For an example of calling this function, see `../unit_tests/replicate/CBIG_MMLDA_brain_behavior_to_doc_wrapper.m`.

----

## Step2: Estimating latent factors with different random initializations
To run MMLDA estimation, you need to specifty number of factors `k`, random seed `r` between`start` and `end` index for different random initializations. The output will be in the folder `${out_dir}/k${k}/r${r}`.

We estimate the following MMLDA model parameters:
* `[iteration].gamma`: regard it as Pr(Factor | Patient), factor compositions
* `[iteration].beta1`: regard it as Pr(Voxel | Factor), probabilistic atrophy maps
* `[iteration].beta2`: regard it as Pr(Score | Factor), probabilistic cognitive deficits
* `[iteration].other`: alpha controlling the distribution of factor compositions  

PS: "regard it as" means they are not exactly the same, you need to some normalization and exponential operation in next step.

To see detailed usage, type `CBIG_MMLDA_est.sh` in terminal. For an example of calling this function, see `../unit_tests/replicate/CBIG_MMLDA_runMMLDA_est_wrapper.sh` or `../examples/CBIG_MMLDA_example_wrapper.m`

----

## Step3: Visualizing factors 
To visualize factors, we will
1. Pick the initialization with maximum log-likelihood
2. Copy final results from MMLDA estimation folder to visualization output folder
3. Plot the climbing of log-likelihood as the best run iterates. The log-likelihood 
   should "go flat" if MMLDA has converged
4. Plot the correlation of each run with best run and the runs are sorted based on
   likelihood. If many initializations yield similar (high correlations) results, then it is enough.
5. Write Pr(Factor | Subject), i.e., normalized gamma in MMLDA, to txt, where each row is a subject.
6. Write Pr(Voxel | Factor), i.e., beta1 in MMLDA estimation, to nifti file.
7. Overlay input volume on underlay volumn (e.g., MNI template) with Hot colormap.
   Then, take screeshots of different slices of the volume and save it as images.
   Please make sure that your VNC resolution is 2000x1000.
8. Visualize behavior deficits of each factor with bar plot. 
9. Project brain from MNI space to fsaverage space and visualize it with Freesurfer.

To see detailed usage, type `help CBIG_visualize_factors` in matlab. For an example of calling this function, see `../unit_tests/replicate/CBIG_MMLDA_visualize_factors_wrapper.sh` or `../examples/CBIG_MMLDA_example_wrapper.m`

----

## Step4: Inferring new subjects with estimated parameters
To run MMLDA inference, you need to specifty number of factors `k`, and estiamted parameters `${visualize_dir}/k${k}/r${r}/final` for the final model. 

We infer the following parameters:
* `[out_name]-gamma.dat`: regard it as Pr(Factor | Patient), factor compositions

To see detailed usage, type `CBIG_MMLDA_inf.sh` in terminal. For an example of calling this function, see `../unit_tests/replicate/CBIG_MMLDA_runMMLDA_inf_wrapper.sh` or `../examples/CBIG_MMLDA_example_wrapper.m`


