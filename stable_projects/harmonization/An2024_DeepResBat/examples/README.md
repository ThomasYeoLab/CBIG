# Examples for DeepResBat harmonization

## References

-   An, L., Zhang, C., Wulan, N., Zhang, S., Chen, P., Ji, F., Ng, KK., Chen, C.,Zhou, J., Yeo, B.T., 2024. [DeepResBat: deep residual batch harmonization accounting for covariate distribution differences](https://doi.org/10.1101/2024.01.18.574145), BioRxiv

---

## Usage

---

These steps shows how to apply DeepResBat harmonization on toy dataset.

### 1. Generate toy data

```
cd $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat

conda activate CBIG_An2024

python -m examples.CBIG_DeepResBat_example --gen_data
```

### 2. Estimate effects of covariates

Our DeepResBat utilizes nonlinear regression trees to estimate effects of covariates. Please refer to section 2.7.1 in manuscript for more details.

```
python -m examples.CBIG_DeepResBat_example --est_covar_effects
```

### 3. Generate covariate-free residuals

The covariate-free residuals are computed by substracting estimated covariates effects from raw data. Please refer to Eq.(12) in manuscript for more details.

```
python -m examples.CBIG_DeepResBat_example --gen_res
```

### 4. HORD search for residual harmonization (optional)

HORD (Ilievski et al., 2017) is a Bayes framework for optimal hyperparameters search. We applied HORD search for cVAE (Section 2.7.2) to get best performance. This step is optional since we found cVAE is relatively stable regarding to hyperparameters, users are recommened to try our emprical hyperparameters firstly. HORD search takes about 20 hours on RTX3090 to finish.

```
python -m examples.CBIG_DeepResBat_example --HORD
```

### 5. Harmonize residuals

Our DeepResBat utilizes cVAE to harmonize covariate-free residuals. Please refer to section 2.7.2 in manuscript for more details.

```
python -m examples.CBIG_DeepResBat_example --harm_res
```

### 6. Add covariates effects back to harmonized residuals

The estimated covariates effects will be added back to harmonized residuals to get harmonized data. Please refer to Eq.(14) in manuscript for more details.

```
python -m examples.CBIG_DeepResBat_example --add_back
```

### 7. Run MANOVA for evaluation

Run MANOVA as harmonization evaulation. Please refer to section 2.9.2 in manuscript for more details.

```
python -m examples.CBIG_DeepResBat_example --manova
python -m examples.CBIG_DeepResBat_example --print_results
```

You should see very close results as follow:

| Beta | -log(P)  | PillarTrace |
| ---- | -------- | ----------- |
| AGE  | 0.173499 | 0.416988    |
| SEX  | 0.433071 | 0.455762    |
| MMSE | 4.039143 | 0.625207    |
| ICV  | 0.940041 | 0.499431    |

---

### Clean up

Run the following commands to clean up after running experiments

```
rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/examples/data

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/examples/checkpoints

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/examples/results

rm -rf $CBIG_CODE_DIR/stable_projects/harmonization/An2024_DeepResBat/examples/__pycache__
```

---

## Bugs and Questions

Please contact Lijun An at anlijuncn@gmail.com and Chen Zhang at chenzhangsutd@gmail.com
