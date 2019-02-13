# Examples

Toy examples of the variance component model and the kernel ridge regression are given here. 
Due to the consideration of file size and the data restriction, here we only utilized 2 subjects (in total 4 sessions) from the CoRR-HNU datasets:
+ **CoRR-HNU:**
  
  `$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0?/subj0?_sess?`
  
Since we cannot do cross-validation with only 2 subjects, each session was treated as a single subject. Then the resting-state functional connectivity matrices of the new 4 "subjects" were duplicated (but with some additional noise) to compose 8 fake "subjects". We then made up two behavioral measures (see `./input/Faked_CSV1.csv`) and two covariates (see `./input/Faked_CSV2.csv`) for the 8 "subjects". Note that this is just a toy example. None of the two models can be applied on only 8 subjects in real studies.

## Variance component model

`./scripts/CBIG_LiGSR_example_Variance_Component.sh` runs though the entire toy example for the variance component model. It first prepares the resting-state functional connectivity matrices, and the functional similarity matrix, and then calls `$CBIG_CODE_DIR/external_packages/LME/Morphometricity_1RandEffect.m` to perform the estimation of explained variances on the full set. Later this function generates 10 delete-2 jackknife samples, and estimates the explained variances in each of these jackknife samples.

Note that this example is only for the goal of showing users about how to use the general variance component model functions. Since there are only 8 faked subjects, the estimated variances are not meaningful at all.

## Kernel ridge regression

`./scripts/CBIG_LiGSR_example_KRR.sh` runs through the entire toy example for kernel ridge regression method. It first prepares the input matrices and parameters required by the kernel ridge regression workflow code, and then calls `$CBIG_CODE_DIR/utilities/matlab/predictive_models/KernelRidgeRegression/CBIG_KRR_workflow.m` to perform the prediction procedure. You can run this script directly on your linux terminal without any parameter passed in.
