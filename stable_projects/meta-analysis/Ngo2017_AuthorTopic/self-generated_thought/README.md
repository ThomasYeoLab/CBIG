This folder contains wrapper functions necessary to reproduce the results on cognitive components of self-generated thought (Fig. 4)

# Initialization
From Matlab, run:
```
initCVB; % add the necessary paths
```

# Inference
To peform inference with the number of cognitive components K ranging from 1 to 4, with 100 random initialization for each K, run in Matlab the following command:
```
CBIG_RunInferenceWrapper(1:4, 1:100);
```
The 100 estimates from the 100 initializations for each value of K would be stored under `outputs/K1` to `outputs/K4`

# Get the best estimate
To get the estimate with the highest variational bound from the 100 estimates for each value of K ranging from 1 to 4, run in Matlab the following command from this directory:
```
CBIG_ComputeBestSolutionWrapper(1:4, 1:100);
```
The best solution for each value of K would be stored under `outputs/best_solution`

# Compute Bayesian Information Criterion (BIC)
To compute the BIC value corresponding to each value of K ranging from 1 to 4, run the followign command in Matlab from this directory:
```
CBIG_ComputeBICWrapper(1:4)
```
This would produce the BIC plot in Fig. 4a.

# Visualize the cognitive components
To project the cognitive components to the brain surface for visualization, run the following commands in Matlab from this directory:
```
K = 2;
CBIG_VisualizeComponentsWrapper(K);
```
This would produce the surface visualization of the 2 cognitve components in Fig. 4b. The images are available under `outputs/best_solution/alpha100_eta0.01/figures`. The brain volumes corresponding to the two cognitive components are available udner `outputs/BIC/ComponentsSmoothness/`

# Get task loading of each component
The probability of each task recruiting a component is available in the `theta` matrix of the struct `params` stored in `outputs/best_solution/K2/alpha100_eta0.01/BestSolution_K002.mat`.
The 7 rows of the `theta` matrix correspond to 7 tasks involving self-generated thought:  `Autobiographical Memory`, `Theory of Mind (story)`, `Theory of Mind (nonstory)`, `Deactivation`, `Navigation`, `Moral Cognition`, `Narrative Comprehension` respectively.

See folder `input_data` for the experiments included under each task.
