# Unit test for Meta-matching (MM)

The unit test of MM contains two parts:
1. Kernel ridge regression (MATLAB)
3. Deep neural network (Python)

For MATLAB code, we have tested in R2018b.

----

## References
+ He, T., An, L., Chen, P., Chen, J., Feng, J., Bzdok, D., Holmes, A.J., Eickhoff, S.B. and Yeo, B.T., 2022. [**Meta-matching as a simple framework to translate phenotypic predictive models from big to small data**](https://doi.org/10.1038/s41593-022-01059-9), Nature Neuroscience 25, 795-804.

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
