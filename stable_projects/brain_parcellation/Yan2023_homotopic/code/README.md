Overview
====

The code utilized in this study includes 2 key steps:

**Step 1: Generate premultiplied fMRI data matrix**

- Directory: `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Yan2023_homotopic/code/step1_generate_fmri_input`
- Wrapper function: `CBIG_hMRF_generate_premultiplied_matrix.m`

**Step 2: Generate hMRF parcellation**
   
- Directory: `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Yan2023_homotopic/code/step2_generate_parcellation`
- Wrapper function: `CBIG_hMRF_wrapper_generate_homotopic_parcellation.m`

**Examples**

We provide detailed examples on how to generate the premultiplied matrix and run the parcellation code under the `examples` folder. **We highly recommended the users to go through the example tutorial first**. Refer to the README.md file under the `examples` folder for more details.

Usage
====

Our parcellation was generated on the `fsaverage6` surface space. The code by so far was only extensively tested on the `fsaverage6` surface space. In the future, we may extend the code into other surface spaces like `fs_LR`.

### Step 1: Generate premultiplied fMRI data matrix
----

We first concatenate the fMRI time courses, normalize and then do a inner product in order to obtain a premultiplied fMRI data matrix. This trick could expediate the processing for the optimization process in step 1. For mathematical details, refer to our supplementary S2.

The description of input arguments are detailed in the function header of the function `CBIG_hMRF_generate_premultiplied_matrix.m`. You can also refer to the input files under the `examples` folder to get a sense of what each of these input files should look like.

To run the wrapper function, start Matlab and cd to the parent folder of the wrapper function, specify the input parameters, and run 
```
CBIG_hMRF_generate_premultiplied_matrix(output_dir, lh_fmri_fullpath_txt, rh_fmri_fullpath_txt, mesh_type, censor_file_fullpath)
```

It may take a few hours to generate the premultiplied matrix, depending on the number of subjects being processed as well as your hardware. On our server, to generate the matrix based on the resting-state data of 1479 GSP subjects takes about 8 hours to run on 1 core. You can always check the output log files from the Matlab command to track the progress of computation.

**Output**

+ `<output_dir>/premultiplied_matrix_single.mat`

and a few other intermediate files, which can be removed if the final matrix has already been generated:

+ `<output_dir>/time_mat/lh_time_matrix.mat`
+ `<output_dir>/time_mat/rh_time_matrix.mat`
+ `<output_dir>/mult_mat/lh_mult_matrix.mat`
+ `<output_dir>/mult_mat/rh_mult_matrix.mat`

**Step 2: Generate hMRF parcellation**
----

The premultiplied matrix is the key input for the parcellation generation algorithm. Other fixed input data such as the gradient map (Gordon2016) and the xyz coordinates of each vertex given the surface mesh (which are precomputed or transformed to the appropriate format) are saved under the `./step2_generate_parcellation/input` folder.

The parcellation algorithm utilizes the Matlab version of the graphcut algorithm (Delong2012). The code can be found here: `https://github.com/ThomasYeoLab/CBIG/tree/master/external_packages/matlab/default_packages/graph_cut`.


## Parameter types and settings
Other than the input data, the parameters are also critical for the final parcellation output. We use the `CBIG_hMRF_set_params` function to output a Matlab struct, which contains input information as well as various parameters. In `CBIG_hMRF_set_params`, there are 3 types of parameters:
- `Non-optional parameter`: those that the user has to specify, including `lhrh_avg_file` and `output_dir`
    + `lhrh_avg_file`: path to the output final whole brain premultiplied matrix in `Step 1`.
    + `output_dir`: desired output path for parcellations generated in `Step 2`.
- `Configurable parameter`: those that are optional and could be configured by the user. 
    + If unspecified, default values would beused.
    + We do not have fixed recommended values, since they are highly variable and dependent on input data dimensionality as well as parcellation resolution. Refer to `CBIG_hMRF_set_params` function header for further details of hyperparameter configuration.
- `Non-configurable parameter`: those that are often fixed. If the user had to configure such a hyperparameter, the user can modify the code in `CBIG_hMRF_set_params`.
Refer to `CBIG_hMRF_set_params` function header for further details.

## Running multiple random initializations
Since the final parcellation result is largely dependent on the random initialization as well, we would recommend the user to run multiple random initializations in parallel and pick the best parcellation among all candidates. To configure this, you should utilize `start_seed` and `num_rand_inits` from `Configurable parameter`. A few examples:
- `start_seed` = 1, `num_rand_inits` = 5: 5 random initializations, based on random seeds {1, 2, ..., 5}.
- `start_seed` = 100, `num_rand_inits` = 5: 5 random initializations, based on random seeds {100, 101, ..., 104}.
- `start_seed` = 100, `num_rand_inits` = 100: 100 random initializations, based on random seeds {100, 101, ..., 199}.

However, since you may want to speed up the generation process and parallelize the parcellation jobs (e.g., one random initialization per job), you can also set `num_rand_inits` to 1, and use a different `start_seed` for each job.

To pick the best result among the seeds, this is what we did in the paper during the benchmarking step: we pick the seed based on lowest `results.E`, which corresponds to the final negative log likelihood of the cost function.

## How to run
Lastly, to run the wrapper function, start Matlab and cd to the parent folder of `CBIG_hMRF_wrapper_generate_homotopic_parcellation.m`, specify the input parameters, and run 
```
CBIG_hMRF_wrapper_generate_homotopic_parcellation(your_premultiplied_matrix_dir, your_output_dir, 'num_clusters', 400, ...)
```

**Output**

- **`<output_dir>`**
    + `<output_dir>/convergence_log` contains the log information of algorithm convergence for each seed.
    + `<output_dir>/results` contains the final parcellation result for all seeds. Each `.mat` file consists of 2 parts:
        - `params`: contains the input and hyperparameter information. Please look up `CBIG_hMRF_set_params.m` for the meaning of each parameter.
        - `results`: contains the final parcellation output, including the final output information from the Graphcut algorithm.
            + `lh_label`: the full cortical label of the left hemisphere, with zeros indicating the medial wall vertices.
            + `rh_label`: the full cortical label of the right hemisphere, with zeros indicating the medial wall vertices.
            + `full_label`: the full cortical label across both hemispheres, with zeros indicating the medial wall vertices.
            + `initial_full_label`: the full cortical label across both hemispheres saved right after parcellation initialization.
            + `D`: final data cost from the cost function.
            + `S`: final smooth cost from the cost function.
            + `E`: final negative log likelihood (data cost + smoothcost).
            + `kappa`: final kappa in the cost function.
            + `initial_d`(only available if params.decrease_d = true): the initial input hyperparameter d.
            + `final_d`(only available if params.decrease_d = true): the final hyperparameter d after decrease-d algorithm.
            + `initial_c`(only available if params.cecrease_c = true): the initial input hyperparameter c.
            + `final_c`(only available if params.cecrease_c = true): the final hyperparameter c after decrease-c algorithm.
    + `<output_dir>/seed_<?>` contains all intermediate parcellation files for each individual seed. Depending on your hyperparameter configuration, you should see at least some the following files:
        - `?_pre_hyperpara_autotune.mat` is the intermediate parcellation before running any of the hyperparameter-autotune algorithm.
        - `?_after_decrease_d.mat` is the intermediate parcellation after running decrease-d algorithm.
        - `?_after_decrease_c.mat` is the intermediate parcellation after running decrease-c algorithm.
        - `?_after_increase_tau.mat` is the intermediate parcellation after running increase-tau algorithm.
        

