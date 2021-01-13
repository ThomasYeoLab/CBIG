# Examples of generating the speed-up version of RSFC gradient map

The examples contain two tasks:
1. Generating RSFC gradient map
   - CoRR-HNU data in fsaverage6 surface space: 1 subjects, 2 runs
2. Estimating diffusion embedding matrices of the RSFC gradient map
   - CoRR-HNU data in fsaverage6 surface space: 1 subjects, 2 runs
   
Note that for simplicity, the example data of task 1 and task 2 are the same. 

----

References
==========
+ Kong et al. Individual-Specific Areal-Level Parcellations Improve Functional Connectivity Prediction of Behavior. Under review.

+ Laumann et al. (2015). Functional System and Areal Organization of a Highly Sampled Individual Human Brain. Neuron 87, 657–670

+ Gordon et al. (2016). Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations. Cerebral Cortex 26, 288–303.

+ Margulies et al. (2016). Situating the default-mode network along a principal gradient of macroscale cortical organization. Proc. Natl. Acad. Sci. U.S.A. 113, 12574–12579.

----

Data
====
In this example, we will perform the algorithms on the resting-state fMRI of CoRR-HNU datasets:
+ **CoRR-HNU:**
  
  `$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0?`

----

Run
====

### Generating input data
----

In the terminal, specify the output directory and call the example data generation script:
```
$CBIG_CODE_DIR/utilities/matlab/speedup_gradients/examples/CBIG_SPGrad_create_example_input_data.sh <output_dir>
```
Note that in your `<output_dir>`, there will be a folder:
+ `<output_dir>/data_list`

We will then provide examples for each task separately.

----

### Generating RSFC gradient map
----

After [**Generating input data**](#generating-input-data), the `data_list` folder under `<output_dir>`:

+ `<output_dir>/data_list/fMRI_list/?h_sub1.txt`

contains the left and right hemisphere fMRI data in `fsaverage6` for subject 1. 

+ `<output_dir>/data_list/censor_list/sub1.txt`

contains the censor lists for subject 1. The outlier frames will be denoted as 0s. 

Start Matlab, in Matlab command window, run the following commands to generate RSFC gradient map for subject 1:

```
output_dir = <output_dir>;
lh_fMRI_files = '<output_dir>/data_list/fMRI_list/lh_sub1.txt';
rh_fMRI_files = '<output_dir>/data_list/fMRI_list/rh_sub1.txt';
censor_files = '<output_dir>/data_list/censor_list/sub1.txt';
lh_ind_surf = 'NONE';
rh_ind_surf = 'NONE';
sub_FC = '100';
sub_verts = '200';
medial_mask = 'NONE';
CBIG_SPGrad_RSFC_gradients(lh_fMRI_fies, rh_fMRI_fies, censor_files, lh_ind_surf, rh_ind_surf, 'fsaverage6', medial_mask, sub_FC, sub_verts, output_dir);
```

The RSFC gradient map of subject 1 is stored under:

+ `<output_dir>/gradients_edge_density.dtseries.nii`

----

### Estimating diffusion embedding matrices of the RSFC gradient map
----

We will first downsample the generated RSFC gradient map from previous step, and compute the geodesic gradient distance matrices:

run the following commands:

```
CBIG_SPGrad_generate_gradient_matrix('fsaverage6', 'NONE', '3.2', output_dir);
```

The estimated gradient geodesic distance matrix will be saved as:
+ `<output_dir>/lh_gradient_distance_matrix.npy`
+ `<output_dir>/rh_gradient_distance_matrix.npy`

We will then apply the diffusion embedding algorithm on the gradient distance matrices.

run the following commands:

```
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
num_component = '100';
cmd = ['sh ' CBIG_CODE_DIR '/utilities/matlab/speedup_gradients/CBIG_SPGrad_diffusion_embedding.sh ' output_dir ' ' num_component];
system(cmd); 
```

The downsampled difffusion embeding results will be saved in output directory: 

+ `<output_dir>/lh_emb_<num_component>_distance_matrix.mat`
+ `<output_dir>/rh_emb_<num_component>_distance_matrix.mat`

Finally, we will upsample the diffusion embedding matrcies back to the orginal resolution.

run the following commands:

```
CBIG_SPGrad_upsample_embed_matrix('fsaverage6', 'NONE', 100, output_dir)
```

The upsampled difffusion embeding results will overwrite the downsampled version in output directory: 
+ `<output_dir>/lh_emb_<num_component>_distance_matrix.mat`
+ `<output_dir>/rh_emb_<num_component>_distance_matrix.mat`

The user can check and run the example wrapper script:

+ `$CBIG_CODE_DIR/utilities/matlab/speedup_gradients/examples/CBIG_SPGrad_example_wrapper.m`

The user can compare their results with the reference results:

+ `$CBIG_CODE_DIR/utilities/matlab/speedup_gradients/examples/ref_results`