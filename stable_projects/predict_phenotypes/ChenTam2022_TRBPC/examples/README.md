# Examples of Task-Rest Behavioral Prediction in Children (TRBPC)
In this example, we faked functional connectivity (predictor), behavioral data (target), as well as age and sex (confounds) for 100 subjects.
Using these data, we can check if our codes can run correctly. Following modules are checked:
* single-kernel kernel ridge regression
* multi-kernel kernel ridge regression
* linear ridge regression
* predictive feature matrices for the above prediction models

# Usage
* to run all the examples sequentially, run `CBIG_TRBPC_example_wrapper.m`
* to check if your results are correct, run `CBIG_TRBPC_check_example_results.m`