# Examples for comparison of Kernel Regression and Deep Neural Network for RSFC based behavioral prediction (KRDNN)

The examples of KRDNN contain three parts:
1. Kernel regression with cross validation (MATLAB)
2. Kernel regression with normal training, validation and testing (MATLAB)
3. Deep neural network (Python)

For MATLAB code, we have tested in R2014a and R2018b.

----

References
==========
+ He, T., Kong, R., Holmes, A., Nguyen, M., Sabuncu, M., Eickhoff, S.B., Bzdok, D., Feng, J. and Yeo, B.T., 2019. [**Deep Neural Networks and Kernel Regression Achieve Comparable Accuracies for Functional Connectivity Prediction of Behavior and Demographics**](https://www.biorxiv.org/content/10.1101/473603v1), under review.

----

Data
====
In this example, the data of part 1 and part 2 are synthetic data. The data of part 3 are [ADHD data from nilearn](https://nilearn.github.io/modules/generated/nilearn.datasets.fetch_adhd.html).

----

Run
====

### Kernel regression with cross validation (MATLAB)
----
Run the following code in matlab from this `example` folder:
```MATLAB
CBIG_KRDNN_KR_CV_example
```
You can check the output result in `output_kr_cv_example/final_result.mat` against the reference result in: `results/KR_CV_matlab_r2014a/final_result.mat` or `results/KR_CV_matlab_r2018b/final_result.mat` (There are little differences between matlab R2014a and R2018b, which is might due to the unstable of very small example dataset).

----

### Kernel regression with normal training, validation and testing (MATLAB)
----
Run the following code in matlab from this `example` folder:
```MATLAB
CBIG_KRDNN_KR_TVT_example
```
You can check the output result in `output_kr_tvt_example/final_result.mat` against the reference result in: `results/KR_TVT/final_result.mat`

----

### Deep neural network (Python)
----
Run the following code in terminal from this `example` folder:
```sh
export PYTHONPATH=${PYTHONPATH}:$PWD/../
python3 CBIG_KRDNN_DNN_example.py
```
You can check the terminal output against the reference output in: `results/DNN_ref_output.txt`. The result for FNN and BrainNetCNN (UK Biobank version) should be the same as the reference result. While the result for GCNN (both version), FNN and BrainNetCNN (HCP version) is going to slightly differ from the reference result due to Keras.

----

Bugs and Questions
====
Please contact He Tong at hetong1115@gmail.com.
