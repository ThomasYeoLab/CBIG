# Unit test for comparison of Kernel Regression and Deep Neural Network for RSFC based behavioral prediction (KRDNN)

The unit test of KRDNN contains two parts:
1. Kernel regression (MATLAB)
3. Deep neural network (Python)

For MATLAB code, we have tested in R2018b.

----

References
==========
+ He, T., Kong, R., Holmes, A., Nguyen, M., Sabuncu, M., Eickhoff, S.B., Bzdok, D., Feng, J. and Yeo, B.T., 2019. [**Deep Neural Networks and Kernel Regression Achieve Comparable Accuracies for Functional Connectivity Prediction of Behavior and Demographics**](https://www.biorxiv.org/content/10.1101/473603v1), under review.

----

Data
====
The data of part 1 are synthetic data. The data of part 2 are [ADHD data from nilearn](https://nilearn.github.io/modules/generated/nilearn.datasets.fetch_adhd.html).

----

Run
====

### Kernel regression
----
Run the following code in matlab:
```MATLAB
runtests('CBIG_KRDNN_unit_test')
```

----

### Deep neural network (Python)
----
Run the following code in terminal:
```sh
export PYTHONPATH=${PYTHONPATH}:$PWD/../
python3 CBIG_KRDNN_unit_test_dnn.py
```

----

Bugs and Questions
====
Please contact He Tong at hetong1115@gmail.com.
