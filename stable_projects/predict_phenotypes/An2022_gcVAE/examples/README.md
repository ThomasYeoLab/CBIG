# Examples of Goal-specific brain MRI harmonization


## References

+ An, L., Chen, J., Chen, P., Zhang, C., He, T., Chen, C., Zhou, J., Yeo, B.T., 2022. [Goal-specific brain MRI harmonization](https://doi.org/10.1016/j.neuroimage.2022.119570), NeuroImage, In press

----

## Usage

### 1. Generate toy data


```
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE

conda activate CBIG_An2022

python -m examples.CBIG_gcVAE_example --stage prepare
```

### 2. Train goal-specific DNN model
```
python -m examples.CBIG_gcVAE_example --model goalDNN --stage train
```

### 3. Train cVAE model
```
python -m examples.CBIG_gcVAE_example --model cVAE --stage train
```


### 4. Train gcVAE model
```
python -m examples.CBIG_gcVAE_example --model gcVAE --stage train
```

### 5. Harmonize data
```
python -m examples.CBIG_gcVAE_example --model cVAE --stage predict

python -m examples.CBIG_gcVAE_example --model gcVAE --stage predict
```


### 6. Evaluate harmonization performance

```
python -m examples.CBIG_gcVAE_example --model goalDNN --stage predict

python -m examples.CBIG_gcVAE_example --model XGBoost --stage train

python -m examples.CBIG_gcVAE_example --stage print
```

You are expected to get similar results as following:

```
Without harmonization:
    MMSE prediction MAE: 6.0008
    Dataset prediction accuracy: 0.9701
-------------------------------------------
cVAE harmonization:
    MMSE prediction MAE: 3.9773
    Dataset prediction accuracy: 0.7612
-------------------------------------------
gcVAE harmonization:
    MMSE prediction MAE: 3.7643
    Dataset prediction accuracy: 0.7761
```

----

### Clean up

Run the following commands to clean up after running experiments

```
rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/examples/data

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/examples/checkpoints

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/examples/results
```

----

## Bugs and Questions
Please contact Lijun An at anlijuncn@gmail.com, Pansheng Chen at chenpansheng@gmail.com and Chen Zhang at chenzhangsutd@gmail.com
