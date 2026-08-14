## Examples

This folder contains a small example that checks whether the analysis code is set
up correctly on your machine.

`CBIG_EIAD_example.m` runs a reduced version of the ADNI analysis in
`replication/script/CBIG_EIAD_ADNI.m` on the shared participant table
(`ADNI_df.csv`, 302 participants):

1. a GLM between the mean cortical E/I ratio and CSF Abeta42;
2. a region-wise GLM of the E/I ratio against CSF Abeta42, and the Spearman
   correlation of the resulting 68 regional slopes with the
   sensorimotor-association (SA) axis;
3. a mediation analysis of the CSF Abeta42 -> memory relationship with the mean
   cortical E/I ratio as the mediator.

The regional mediation and the surface figures of the full analysis are skipped
to keep the run time short.

### A note on the shared data

To comply with the ADNI data-use agreement, the only real values released in
`ADNI_df.csv` are the **subject IDs** and the **68 harmonized regional E/I
ratios** (`HarmonizedEI1`–`HarmonizedEI68`). The `AB42_CSF` and `PHCMemory`
columns are **synthetic** and are provided only so that this example runs 
end to end.

Because the CSF Abeta42 and memory columns are synthetic, and because the
covariates used in the paper (age, sex, in-scanner motion, and years of
education) are not included, the numbers this example prints are **illustrative
only** and are not the values reported in the paper. To reproduce the published
results, obtain the real `AB42_CSF`, `PHCMemory`, and covariates from ADNI using
the released subject IDs (see the replication instructions).

### Requirements

MATLAB with the Statistics and Machine Learning Toolbox (`fitglm`, `zscore`,
`randsample`).

### Usage

From the `script` folder, run

```matlab
CBIG_EIAD_example
```

in MATLAB. The results are written to `../output/example_output.mat`.

The results are then compared against the reference 
output (`../reference_output/reference_output.mat`) field by
field with a tolerance of 1e-6.

### Expected result

The script prints the GLM coefficient table and the
mediation paths and ends with

```
Passed! Results are the same as the reference outputs.
```

If a field differs by more than the tolerance, the script reports which field it
was; the most likely cause is a difference in MATLAB version or hardware rather
than an error in the setup. Typical run time is under 10 seconds.
