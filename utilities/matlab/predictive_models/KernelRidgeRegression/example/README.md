# Kernel Ridge Regression (KRR) Example



This fold provides a simple and concrete example to let users familiarize our KRR workflow. The original dataset used for this example is a publicly-available air quality dataset: [Air Quality Dataset](https://archive.ics.uci.edu/ml/datasets/air+quality#). In short, our goal is to predict the hourly-averaged concentration of carbon monoxide (CO) in the air using KRR. Note that this KRR example is for the purpose of illustration of the KRR script usage only, please refrain from making any conclusion regarding the air quality dataset per se. This example is performed using 5-fold cross validation with 3 different splits. To shorten the computational time, only the first 15 days of air quality data are used for this example, the air quality dataset can be found at `/input/AirQuality.csv`.  



## How to use the scripts

There are two versions of the scripts available for this example:

1. Setup file preparation option:

   If the user chooses this option, a setup file for KRR will be prepared and generated before passing it to the KRR workflow. To run this version of example scripts, you can execute the following command line into a Linux terminal:

   `$CBIG_CODE_DIR/utilities/matlab/predictive_models/KernelRidgeRegression/example/script/CBIG_KRR_example_wrapper_setup_file.sh <output_dir>`

   You need to specify the path of your desired output directory.

2. Input arguments option:

   If the user chooses this option, all the necessary input arguments will be generated inside the script. To run this version of example scripts, you can run the following command in Matlab:

   `CBIG_KRR_example_wrapper_input_args(output_dir)`

   You need to specify the path of your desired output directory.



## KRR example setup file preparation/input arguments

Prior to running KRR, certain parameters must be defined. For the purpose of this example, the following parameters are specified  as follows:

* `outdir`

  This parameter specifies the path of output directory for each split.

* `sub_fold`

  This parameter contains the information of how the dataset is split. In this example, we treat each row of `AirQuality.csv` as an individual subject.  We ignore certain rows that contain missing data. As a result, there are 312 data points (subjects) in total. The data identifier (data-id) can be found at `/input/data-id.txt`.  

* `y`

  This parameter contains the ground truth value for prediction. In this example, the only prediction is about the hourly-averaged concentration of CO. The values can be found at `/input/CO_groundtruth.csv`. The unit is in mg/m^3.

* `covariates`

  In this example, we treat `date` as a covariate and regress out `date` from the target variable (hourly-averaged concentration of CO). The `date` information for each data point can be found at `/input/date.csv`. The value for `date` range from 1 to 15, indicating that only 15 days of the data are used for this example.

* `feature_mat`

  `feature_mat` is a F x N matrix, where F is the number of features and N is the number of data points. The features used for this example include:

  * PT08.S1(CO): PT08.S1 (tin oxide) hourly averaged sensor response (nominally CO targeted)
  * C6H6(GT): True hourly averaged Benzene concentration in microg/m^3
  * PT08.S2(NMHC): True hourly averaged overall Non Metanic HydroCarbons concentration in microg/m^3 (reference analyzer)
  * NOx(GT): True hourly averaged NOx concentration in ppb (reference analyzer)
  * PT08.S3(NOx): PT08.S3 (tungsten oxide) hourly averaged sensor response (nominally NOx targeted)
  * NO2(GT): True hourly averaged NO2 concentration in microg/m^3 (reference analyzer)
  * PT08.S4(NO2): PT08.S4 (tungsten oxide) hourly averaged sensor response (nominally NO2 targeted)
  * PT08.S5(O3): PT08.S5 (indium oxide) hourly averaged sensor response (nominally O3 targeted)
  * T: Temperature in °C
  * RH: Relative Humidity (%)
  * AH: Absolute Humidity

  Therefore, `feature_mat` for this example has a size of 11 x 312. `feature_mat` can be found at `/input/feature.mat`

* `num_inner_folds`

  The number of inner-loop cross-validation folds is set to 5 in this example.

* `outstem`

  This parameter specifies the output file name for each split. 

* `with_bias`

  This example set this parameter to 1 to include a bias term.

* `ker_param`

  The type of kernel used is the Pearson's correlation.* 

* `lambda_set`

  The regularization parameter lambda is optimized through a grid search ranging from `[0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20]`(default range).

* `metric`

  The tuning metric used for this example is `predictive_COD`.

  

## KRR example reference output

The reference output is located at `$CBIG_CODE_DIR/utilities/matlab/predictive_models/KernelRidgeRegression/example/reference_output`. It contains the KRR results with 5-fold cross validation with 3 different splits. Each directory (1, 2, 3) contains the output from each split. The final output is saved as `final_result_CO_[split_seed].mat.`

By running the script (setup file option or input argument option), it automatically compares the generated output with the reference output. An error message will be displayed if the optimal accuracies differ too much from the reference accuracies.



## Bugs and Questions

Please contact Zhang Shaoshi at shaoshi.zhang@u.nus.edu
