## Unit test of E/I ratio and AD pathophysiology (EIAD)

The unit test of EIAD runs the example provided in the `examples` folder, which
performs:

1. a GLM between the mean cortical E/I ratio and CSF Abeta42;
2. a region-wise GLM of the E/I ratio against CSF Abeta42, and the Spearman
   correlation of the resulting 68 regional slopes with the
   sensorimotor-association (SA) axis;
3. a mediation analysis of the CSF Abeta42 -> memory relationship with the mean
   cortical E/I ratio as the mediator.

The example writes its results to `examples/output/example_output.mat` and
compares them field by field against `examples/reference_output/reference_output.mat`
with a tolerance of 1e-6. The test passes if every field agrees. The output
folder is removed afterwards.

If `replace_unittest_flag` is set to 1, the reference output is replaced by the
results of this run instead of being compared against.

See README.md under the `examples` folder for more details. Notice that the unit
test is **for CBIG lab only**.

## Usage

```matlab
runtests('CBIG_EIAD_unit_test')
```
