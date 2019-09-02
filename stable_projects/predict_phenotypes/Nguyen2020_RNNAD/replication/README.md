# Replicating experiment predicting Alzheimer's Disease progression using recurrent neural network (RNN) 

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

- After setting up the Python environment, run the *CBIG_RNN_replication.sh* file
- The result should be as followed
    - RNN on test set, first fold:
    - mAUC 0.9617 bca 0.9141 adasMAE 3.8233 ventsMAE 0.002080
    - RNN on test set, second fold:
    - mAUC 0.9574 bca 0.8912 adasMAE 4.2712 ventsMAE 0.001265
