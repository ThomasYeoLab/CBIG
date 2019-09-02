# Example of predicting Alzheimer's Disease progression using recurrent neural network (RNN) 

References
====

+ Nguyen, M., Sun, N., Alexander D.C., Feng J., Yeo B.T.T., 2018. [**Modeling Alzheimer’s disease progression using deep recurrent neural networks**](https://doi.org/10.1109/prni.2018.8423955), PRNI, 2018.
+ Nguyen, M., He T., An L., Alexander D.C., Feng J., Yeo B.T.T., 2019. **Predicting Alzheimer’s disease progression
using deep recurrent neural networks**, under review.

----

Data
====

- Running this example requires the data from the [Alzheimer's Disease Neuroimaging Initiative](http://adni.loni.usc.edu) (ADNI)
- After getting the spreadsheets from ADNI, follow the [instructions](https://github.com/noxtoby/TADPOLE/blob/master/TADPOLE_readme.txt) to generate the *TADPOLE_D1_D2.csv* file
- This file is the standard dataset of the [TADPOLE Challenge 2017](http://tadpole.grand-challenge.org)
- Copy the *TADPOLE_D1_D2.csv* file generated to the *data* folder


----

Run
====

- After setting up the Python environment, run the *CBIG_RNN_example.sh* file
- The result should be as followed
    - RNN on validation set:
    - mAUC 0.9463 bca 0.8930 adasMAE 3.9777 ventsMAE 0.002248
    - RNN on test set:
    - mAUC 0.9335 bca 0.8772 adasMAE 4.2398 ventsMAE 0.002389
    - Constant prediction baseline on test set
    - mAUC 0.8799 bca 0.8772 adasMAE 4.8312 ventsMAE 0.002712
    - SVM baseline on test set
    - mAUC 0.9182 bca 0.8307 adasMAE 4.9428 ventsMAE 0.002334
