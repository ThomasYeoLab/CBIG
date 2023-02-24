# Predictive-feature matrix (PFM)

This folder contains functions to compute the predictive network features of predictive models by computing the covariance between the features and the predicted behavior (Haufe et al. 2014).

For predictive models of single-kernel ridge regression, multi-kernel ridge regression, and linear ridge regression using fitrlinear from CBIG `$CBIG_CODE_DIR/utilities/matlab/predictive_models`, corresponding PFM calculation functions are `CBIG_compute_singleKRR_PFM.m`, `CBIG_compute_multiKRR_PFM.m`, `CBIG_compute_LRR_fitrlinear_PFM.m`.
For other predictive models, corresponding PFM calculation function is `CBIG_compute_PFM_general.m`.

----

## Usage

To run `CBIG_compute_singleKRR_PFM.m` (similar to `CBIG_compute_multiKRR_PFM.m` and `CBIG_compute_LRR_fitrlinear_PFM.m`), followings are used as an input.

### Inputs:

   - `feature_file`
     The feature file used in the kernel regression (dim: #features x #subjects).
     The variable storing the #features x #subjects matrix can be any name.

   - `singleKRR_dir/multiKRR_dir/LRR_dir`
     The directory where the single-kernel ridge regression/multi-kernel ridge regression/linear ridge regression results are stored. This is the same directory you used as input argument `outdir` when you perform the kernel regression.

   - `sub_fold_file`
     Full path of the cross-validation data split file. A structure `sub_fold` is assumed to be saved in this file. sub_fold(i).fold_index is a #subjects x 1 binary vector. 1 refers to the corresponding subject is a test subject in the i-th test fold. 0 refers to the corresponding subject is a training subject in the i-th test fold.

   - `score_ind`
     A scalar. The index of the score you want to compute the PFM. Range from 1 to # Target Variables in your single-kernel ridge regression. 

   - `outstem`
     A string that was used to describe the .mat file of the target variable y after regression.

   - `outdir`
     Output directory for the predictive-feature matrices.

### Saved results

When `CBIG_compute_singleKRR_PFM.m` (similar to `CBIG_compute_multiKRR_PFM.m` and `CBIG_compute_LRR_fitrlinear_PFM.m`) is run, one # features by # folds predictive-feature matrix will be saved to the output directory for each behavior score.


To run `CBIG_compute_PFM_general.m`, followings are used as an input.

### Inputs:

   - `feature_file`
     The feature file used in predictive model (dim: #features x #subjects).
     The variable storing the #features x #subjects matrix can be any name.

   - `y_pred_train_file`
     The predicted y for training subjects from the predictive model (dim: #subjects x 1).
     The variable storing the #subjects x 1 matrix can be any name.

### Outputs

   - `PFM`
     The predictive-feature matrix (dim: #features x 1).
     
----
## Updates

- Release v0.27.0 (24/02/2023): Initial release of PFM utilities

----
## Bugs and questions
Please contact Naren Wulan wulannarenzhao@gmail.com.
