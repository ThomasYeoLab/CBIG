# Examples of Meta-matching T1


## References

+ To be added

----

## Usage

The experiments contains 8 steps: \
    1. Generate toy data \
    2. Train an ElasticNet model using example test set (which is called meta-test in our paper) \
    3. Train a 3D CNN Model using example train set (which is called meta-training in our paper) \
    4. Run the trained 3D CNN model using example test set to get outputs\
    5. Train a classical transfer model using example test set \
    6. Train a meta-matching finetune model using example test set \
    7. Train a meta-matching stacking model using example test set \
    8. Print results for 4 methods

Run the following command
```

cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1

source activate Naren2024_MMT1

python -m examples.CBIG_MMT1_example

```


You are expected to get similar results as following:

```
>>>>>>> Elasticnet <<<<<<<<
average correlation: 0.9925 0.9999 0.9999 0.9999 0.9999
average correlation for K = 10        20        50        100       200
average COD: 0.9856 0.9999 0.9999 0.9999 0.9999
average COD for K = 10        20        50        100       200

>>>>>>> Directly transfer <<<<<<<<
average correlation: 0.2982 0.1176 0.1674 0.1081 0.3791
average correlation for K = 10        20        50        100       200
average COD: 0.0685 -0.5483 0.1092 -0.4242 -0.0102
average COD for K = 10        20        50        100       200

>>>>>>> MM finetune <<<<<<<<
average correlation: 0.4367 0.4306 0.4213 0.4200 0.6351
average correlation for K = 10        20        50        100       200
average COD: 0.2268 0.2121 0.0927 0.0829 0.0316
average COD for K = 10        20        50        100       200

>>>>>>> MM stacking <<<<<<<<
average correlation: 0.4426 0.4384 0.43380.4328 0.6255
average correlation for K = 10        20        50        100       200
average COD: 0.1852 0.2021 0.2209 0.2188 0.3248
average COD for K = 10        20        50        100       200

```
Please note, this is just a generated data, please refer to the 
replication folder and the paper for the performance on real dataset.

----

### Clean up

Run the following commands to clean up after running experiments

```
rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/examples/data

rm -rf $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/examples/results
```

----

## Bugs and Questions
Please contact Naren Wulan at wulannarenzhao@gmail.com, Chen Zhang at chenzhangsutd@gmail.com and 
Lijun An at anlijuncn@gmail.com
