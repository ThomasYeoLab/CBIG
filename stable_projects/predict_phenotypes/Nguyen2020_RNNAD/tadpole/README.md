# Replicating CBIG submitted entry to the TADPOLE grand challenge

References
====

+ Nguyen, M., Sun, N., Alexander D.C., Feng J., Yeo B.T.T., 2018. [**Modeling Alzheimer’s disease progression using deep recurrent neural networks**](https://doi.org/10.1109/prni.2018.8423955), PRNI, 2018.
+ Nguyen, M., He T., An L., Alexander D.C., Feng J., Yeo B.T.T., 2020. [**Predicting Alzheimer’s disease progression using deep recurrent neural networks**](https://doi.org/10.1016/j.neuroimage.2020.117203), NeuroImage, 117203.

----

Data
====

- Replicating this requires the data from the [Alzheimer's Disease Neuroimaging Initiative](http://adni.loni.usc.edu) (ADNI)
- After getting the spreadsheets from ADNI, follow the [instructions](https://github.com/noxtoby/TADPOLE/blob/master/TADPOLE_readme.txt) to generate the *TADPOLE_D1_D2.csv* file
- This file is the standard dataset of the [TADPOLE Challenge 2017](http://tadpole.grand-challenge.org)
- Copy the *TADPOLE_D1_D2.csv* file generated to the *data* folder


----

Run
====

- Set up the Python environment using *Anaconda* and the *replication/config/CBIG_RNN_python_env.yml* file
- After setting up the Python environment, run the *CBIG_RNN_ensemble_prediction.sh* file
- First, the script generates the necessary data to make prediction
- Second, it loads 4 sets of pre-trained weights and makes 4 predictions for the same set of test subjects.
- The pre-trained weights file (*model_weights.pt*) and the model hyperparameters file (*config.json*) are in the *save.seed0*, *save.seed1*, *save.seed2*, *save.seed3* folders.
- Finally, it combines the 4 predictions by averaging.
- The final output file is named *TADPOLE_Submission_D4Live.csv*.
- Due to difference in CPU and GPU used, this output will not be exactly the same as our submitted entry.
However, it should be very close.
