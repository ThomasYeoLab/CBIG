# Unit test for Meta-matching (MM)

The unit test of MM contains two parts:
1. Kernel ridge regression (MATLAB)
3. Deep neural network (Python)

For MATLAB code, we have tested in R2018b.

----

## References
+ He, T., An, L., Feng, J., Bzdok, D., Eickhoff, S.B. and Yeo, B.T., 2020. [**Meta-matching: a simple approach to leverage large-scale brain imaging datasets to boost prediction of non-imaging phenotypes in small datasets**](https://doi.org/10.1101/2020.08.10.245373), under review.

----

## Data
The data are synthetic data.

----

## Run

### Kernel ridge regression
----
Run the following code in matlab:
```MATLAB
runtests('CBIG_MM_unit_test')
```

----

### Deep neural network (Python)
----
Run the following code in terminal (need gpu to run):
```sh
export PYTHONPATH=${PYTHONPATH}:$PWD/../
python3 CBIG_MM_unit_test_DNN.py
```

----

## Bugs and Questions
Please contact Tong He at hetong1115@gmail.com, Lijun An at anlijun.cn@gmail.com and Pansheng Chen at chenpansheng@gmail.com.
