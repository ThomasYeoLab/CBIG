# Examples of Areal-level Gradient-based Multi-session hierarchical Bayesian model (gMS-HBM)

The examples of Areal-level MS-HBM contain several tasks:
1. Generating diffusion embedding matrices of RSFC gradients
2. Generating profiles, initialization parameters, and spatial radius masks
3. Group priors estimation
4. Individual-level parcellation generation
   
Note that for simplicity, we use the same example data (two CoRR-HNU subjects) for all tasks. However, the initialization parameters estimated in task 2 **won't** be used in task 3, the group priors estimated in task 3 **won't** be used in task 4. We will use the initialization parameters and group priors pre-generated by HCP dataset (40 subjects, 4 sessions). 

----

References
==========
+ Kong R, Yang Q, Gordon E, et al. [Individual-Specific Areal-Level Parcellations Improve Functional Connectivity Prediction of Behavior](https://doi.org/10.1093/cercor/bhab101). Cerebral Cortex. In press.

----

Data
====
In this example, we will perform the algorithms on fsaverage6 surface data of CoRR-HNU datasets:
+ **CoRR-HNU:**
  
  `$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0?/subj0?_sess?`

----

Run
====

### [IMPORTANT NEW FEATURE] A simplified way to generate individual parcellation for a single subject
----

We now provide a simple way to generate individual parcellation for a single subject. The user only needs to specify paths to the fMRI data and censor files, the algorithm will automatically generate the profiles and use group priors to generate individual parcellation for a single subject.

In the terminal, specify the paths:
```
params.project_dir = out_dir;

% user can set their own group prior
%params.group_prior = '<your_own>/priors/Params_Final.mat';

params.censor_list = {fullfile(CBIG_CODE_DIR,'/data/example_data/CoRR_HNU/subj01/subj01_sess1/qc/subj01_sess1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt')};
params.lh_fMRI_list = {fullfile(CBIG_CODE_DIR,'/data/example_data/CoRR_HNU/subj01/subj01_sess1/surf/lh.subj01_sess1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz')};
params.rh_fMRI_list = {fullfile(CBIG_CODE_DIR,'/data/example_data/CoRR_HNU/subj01/subj01_sess1/surf/rh.subj01_sess1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz')};
params.target_mesh = 'fsaverage6';
params.model = 'gMSHBM';
params.w = '50';
params.c = '50';
[lh_labels, rh_labels] = CBIG_ArealMSHBM_parcellation_single_subject(params);
```
We have also provided a wrapper script the user can check:
+ $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/CBIG_ArealMSHBM_example_single_subject.m

The following section is the original example tutorial with more detailed steps. Please note the new feature only can generate parcellation for a single subject. It does not contain the training and validation steps.

### **Generating input data**
----

In the terminal, specify the output directory and call the example data generation script:
```
$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/CBIG_ArealMSHBM_create_example_input_data.sh <output_dir>
```
Note that in your `<output_dir>`, there will be four folders corresponding to four different tasks:
+ `<output_dir>/generate_gradients`
+ `<output_dir>/generate_profiles_and_ini_params`
+ `<output_dir>/estimate_group_priors`
+ `<output_dir>/generate_individual_parcellations`

We will then provide examples for each task separately.

----

### **Example1: generating diffusion embedding matrices of RSFC gradients**
----

This step is mainly needed for gradient-based MSHBM (gMSHBM). If you are only interested in dMSHBM or cMSHBM, you can skip this example.

After [**Generating input data**](#generating-input-data):

+ `<output_dir>/generate_gradients/data_list/fMRI_list/?h_sub?_sess?.txt`

contains the left and right hemisphere fMRI data in `fsaverage6` for subject 1 session 1. Note that our example data only has 1 run for each session, if there are multiple runs, then the corresponding `?h_sub?_sess?.txt` should have multiple columns, where each column corresponds to each run, the columns should be separated by white space.

+ `<output_dir>/generate_gradients/data_list/censor_list/sub?_sess?.txt`

contains the censor lists for subject 1 session 1. The outlier frames will be denoted as 0s. Again, if there are multiple runs, `sub?_sess?.txt` should have multiple columns.

Note that the example data does not have individual surface template files. Therefore, there is no `<output_dir>/generate_gradients/data_list/surface_list`. However, for dataset such as HCP, they do provide the individual surface templates files (e.g. `100206.L.midthickness.32k_fs_LR.surf.gii`). In this case, the user should have a `surface_list` folder. Check the README.md file under `step0_generate_gradient_prior` for more details.

#### 1.1 Generating gradients

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step0_generate_gradient_prior
```

Start Matlab, in Matlab command window, run the following commands to generate gradients for subject 1 session 1, the fMRI data are in `fsaverage6`.

```
project_dir = '<output_dir>/generate_gradients';
CBIG_ArealMSHBM_generate_gradient('fsaverage6',project_dir,'1','1');
```
The diffusion embedding matrices of gradients will be saved into `gradients` folder:

+ `<output_dir>/generate_gradients/gradients/sub1/?h_emb_100_distance_matrix.mat`

----

### **Example2: generating profiles, initialization parameters, and spatial radius masks**
----

After [**Generating input data**](#generating-input-data):

+ `<output_dir>/generate_profiles_and_ini_params/data_list/fMRI_list/?h_sub?_sess?.txt`

contains the left and right hemisphere fMRI data in `fsaverage6` for each subject and each session. Note that our example data only has 1 run for each session, if there are multiple runs, then the corresponding `?h_sub?_sess?.txt` should have multiple columns, where each column corresponds to each run, the columns should be separated by white space.

+ `<output_dir>/generate_profiles_and_ini_params/data_list/censor_list/sub?_sess?.txt`

contains the censor lists for each subject and each session. The outlier frames will be denoted as 0s. Again, if there are multiple runs, `sub?_sess?.txt` should have multiple columns. Furthermore, if the user don't want to mask out the outliers in computing profiles, just leave the `censor_list` folder empty.

#### 2.1 Generating profiles

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step1_generate_profiles_and_ini_params
```

Start Matlab, in Matlab command window, run the following commands to generate profiles for each session of the 2 subjects, the fMRI data are in `fsaverage6`, the ROIs are defined by `fsaverage3`.

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
for sub = 1:2
    for sess = 1:2
        CBIG_ArealMSHBM_generate_profiles('fsaverage3','fsaverage6',project_dir,num2str(sub),num2str(sess),'0');
    end
end
```
The profiles will be saved into `profiles` folder:

+ `<output_dir>/generate_profiles_and_ini_params/profiles/sub?/sess?/?h.sub?_sess?_fsaverage6_roifsaverage3.surf2surf_profile.nii.gz`

If the user only has 1 session and would like to apply it on Areal-level MSHBM, the user need to split the session into multiple sub-sessions to make Areal-level MSHBM work. In practice, to estimate better individual-level areal-level parcellation, we also find it's better to split a scan into multiple sub-sessions if the subject has less than 4 runs of data. We suggest each sub-session have ~5 min data.

In this example, we split session 1 of CoRR-HNU subject 1 into 2 sub-sessions. To do that, run the following commands to generate 2 profiles for each sub-session:

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
CBIG_ArealMSHBM_generate_profiles('fsaverage3','fsaverage6',project_dir,'1','1','1');
```

The profiles will be saved into `profiles` folder:

+ `<output_dir>/generate_profiles_and_ini_params/profiles/sub?/sess?/?h.sub?_sess?_fsaverage6_roifsaverage3.surf2surf_profile_1.nii.gz`

+ `<output_dir>/generate_profiles_and_ini_params/profiles/sub?/sess?/?h.sub?_sess?_fsaverage6_roifsaverage3.surf2surf_profile_2.nii.gz`

#### 2.2 Generating initialization parameters

To generate the initialization parameters, we need to obtain group averaged profiles of the training set. In this example, we use 2 CoRR-HNU subjects with 2 sessions as the training set. Note that the profiles of split sub-sessions will be ignored during averaging. In practice, if the subject has different number of sessions, please set the `num_sess` to be the maximum number of sessions.

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
num_sub = '2';
num_sess = '2';
CBIG_ArealMSHBM_avg_profiles('fsaverage3','fsaverage6',project_dir,num_sub,num_sess);
```

The results will be saved into `profiles/avg_profile` folder:

+ `<output_dir>/generate_profiles_and_ini_params/profiles/avg_profile/?h_fsaverage6_roifsaverage3_avg_profile.nii.gz`

Next, we will use Schaefer 400-area group-level parcellation to generate the initialization parameters with the group average profiles:

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
num_parcels = '400';
CBIG_ArealMSHBM_generate_ini_params_Schaefer('fsaverage3','fsaverage6',num_parcels, project_dir);
```

The results will be saved into `group` folder:

+ `<output_dir>/generate_profiles_and_ini_params/group/group.mat`

contains the estimated initialization parameter `mtc`, which is the group-level functional connectivity profile of parcels. `group.mat` also contains the group-level parcellation `lh_labels` and `rh_labels`.

Please also note that we **will not** use the results generated in this step in the following sections. Instead, we will use a pre-computed initialization parameters from 40 HCP subjects for the next steps. 

#### 2.3 Generating spatial radius mask

In the areal-level MS-HBM model, we introduce spatial radius mask so that the individual-level parcels should be within the radius mask of group-level parcels. In our paper, we use 30mm radius mask.

In this example, we show how to generate radius mask for Schaefer 400-area parcellation in `fsaverage5` for simplicity. Run the following commands to generate the spatial mask in this example:

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
num_parcels = '400';
CBIG_ArealMSHBM_generate_radius_mask_Schaefer(num_parcels, 'fsaverage5', project_dir);
```

The results will be saved into `spatial_mask` folder:

+ `<output_dir>/generate_profiles_and_ini_params/spatial_mask/spatial_mask_fsaverage5.mat`

Note that in this example we don't do it in `fsaverage6` is because it's very slow (~40min). However, we released the pre-computed spatial radius masks for Schaefer parcellations in `fsaverage6` and `fs_LR_32k` here:

+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/spatial_mask`

----

### **Example3: Group priors estimation**
----

Note that we will use pre-computed functional connectivity profiles of 2 CoRR-HNU subject with 2 sessions. There is no need to re-generate the profiles. We will also use pre-computed initialization parameters and spatial radius masks for Schaefer 400-area parcellation.

**Important Note:**

In practice, we use ~40 subjects as training set to estimate group priors. To speed things up, our group prior estimation scripts were written as parent-child relationship. This means that if there are `S` training subjects, the user need to run 1 parent script and `S` children scripts in parallel. If the user would like to train their own group priors with their own data, they need to have a multi-core server which can submit `1+S` jobs which can run in parallel.

In this example, we will estimate group priors with 2 subjects and stop the estimation procedure at the 2nd iteration to save time. In practice, we suggest user to specify the maximum iteration as 10 or leave it empty, the script normally will stop when it converges.

Please also note that we have released our pre-estimated group priors for Schaefer parcellations in `fsaverage6` and `fs_LR_32k` spaces. We suggest the user directly use our pre-estimated group priors to avoid high computational cost in the training step. These group-priors were estimated by 40 HCP subjects. We have tested them in Midnight Scnanning Club (MSC) dataset for `fs_LR_32k` and Adolescent Brain Cognitive Development (ABCD) dataset for `fsaverage6`.

After [**Generating input data**](#generating-input-data),

the `group` folder under `<output_dir>/estimate_group_priors`:
+ `<output_dir>/estimate_group_priors/group/group.mat`

contains the **pre-computed** initialization parameters for Schaefer 400-area parcellation.

The `spatial_mask` folder under `<output_dir>/estimate_group_priors`:
+ `<output_dir>/estimate_group_priors/spatial_mask/spatial_mask_fsaverage6.mat`

contains the **pre-computed** 30mm spatial radius mask for Schaefer 400-area parcellation.

The `profile_list/training_set` folder under `<output_dir>/estimate_group_priors`:
+ `<output_dir>/estimate_group_priors/profile_list/training_set/lh_sess?.txt`
+ `<output_dir>/estimate_group_priors/profile_list/training_set/rh_sess?.txt`

contain functional connectivity profiles of 2 CoRR-HNU subjects with 2 sessions. These files are also **pre-generated**.

The `gradient_list/training_set` folder under `<output_dir>/estimate_group_priors`:
+ `<output_dir>/estimate_group_priors/gradient_list/training_set/gradient_list_lh.txt`
+ `<output_dir>/estimate_group_priors/gradient_list/training_set/gradient_list_rh.txt`

contain diffusion embedding matrices of RSFC gradients of 2 CoRR-HNU subjects generated by 2 sessions. These files are also **pre-generated**.

There is a hyperparameter `beta` for gMSHBM and cMSHBM, which is the weight for the gradient-based/xyz coordinates prior. Larger `beta` indicates stronger weight for the gradient-based/xyz coordinates prior. In practice, the user needs to train model with different `beta` and use validation set to find the optimal `beta`. In this example, we set `beta` to be `5` and we will train with gMSHBM.

The user need to open 3 Matlab sessions to run 1 parent script and 2 chilren scripts. For each Matlab session:

```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step2_estimate_priors
```

In the first Matlab session, run the parent script:
```
project_dir = '<output_dir>/estimate_group_priors';
tmp_dir = '<output_dir>/estimate_group_priors/tmp_results';
beta = '5';
CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_parent(project_dir,'fsaverage6','2','2','400',beta,tmp_dir,'2');
```

In the second Matlab session, run the child script for subject 1:
```
project_dir = '<output_dir>/estimate_group_priors';
tmp_dir = '<output_dir>/estimate_group_priors/tmp_results';
CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_child('1', 'fsaverage6', tmp_dir);
```

In the third Matlab session, run the child script for subject 2:
```
project_dir = '<output_dir>/estimate_group_priors';
tmp_dir = '<output_dir>/estimate_group_priors/tmp_results';
CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_child('2', 'fsaverage6', tmp_dir);
```

The results of each iteration will be saved into `priors` folder, the results of each iteration will be saved as `Params_iteration?.mat`, which contains a struct variable `Params`. The final estimated group priors should be saved as `Params_Final.mat`.

+ `<output_dir>/estimate_group_priors/priors`

The estimated group priors include:
1) Inter-subject functional connectivity variability -- `Params.epsil`
2) Group-level connectivity profiles for each network -- `Params.mu`
3) Intra-subject functional connectivity variability -- `Params.sigma`
4) Spatial prior which denotes the probability of each network occurring at each location -- `Params.theta`

The user can compare their estimated `Params_iteration?.mat` with 
+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_MSHBM/examples/ref_results/estimate_group_priors/priors/gMSHBM/beta5/Params_iteration?.mat`

----

### **Example4: Individual-level parcellations generation**
----

Note that the group priors estimated in [**Group priors estimation**](#Example3:-Group-priors-estimation) won't be used to generate individual parcellations. We will use the group priors pre-computed by using HCP dataset (40 subjects, 4 sessions), which are saved in:

+ `<output_dir>/generate_individual_parcellations/priors/gMSHBM/beta5/Params_Final.mat`


After [**Generating input data**](#generating-input-data),

the `priors` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/priors/gMSHBM/beta5/Params_Final.mat`

contains the **pre-computed** group priors for Schaefer 400-area parcellation with `beta = 5`.

The `spatial_mask` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/spatial_mask/spatial_mask_fsaverage6.mat`

contains the **pre-computed** 30mm spatial radius mask for Schaefer 400-area parcellation.

The `profile_list/validation_set` and `profile_list/test_set` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/profile_list/validation_set/lh_sess?.txt`
+ `<output_dir>/generate_individual_parcellations/profile_list/validation_set/rh_sess?.txt`
+ `<output_dir>/generate_individual_parcellations/profile_list/test_set/lh_sess?.txt`
+ `<output_dir>/generate_individual_parcellations/profile_list/test_set/rh_sess?.txt`

contain functional connectivity profiles of 2 CoRR-HNU subjects with 4 sub-sessions. These files are also **pre-generated**.

The `gradient_list/validation_set` and `gradient_list/test_set` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/gradient_list/validation_set/gradient_list_lh.txt`
+ `<output_dir>/generate_individual_parcellations/gradient_list/validation_set/gradient_list_rh.txt`
+ `<output_dir>/generate_individual_parcellations/gradient_list/test_set/gradient_list_lh.txt`
+ `<output_dir>/generate_individual_parcellations/gradient_list/test_set/gradient_list_rh.txt`

The `data_list/validation_set` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/data_list/validation_set/fMRI_list/?h_sub1.txt`

contain fMRI data of CoRR-HNU subject 1. These fMRI files will be used to evaludate parcellation quality during validation.

#### 4.1 Validation set

We will still use gMSHBM as an example. This algorithm contains three parameters:

+ `beta`: The weight of gradient-based prior.
+ `w`: The weight of group spatial prior `Params.theta`.
+ `c`: The weight of MRF smoothness prior.

These three parameters can be selected by validation set. Assuming subjects in validation set has `T` sessions, we normally use `T1` sessions to generate the individual parcellations for each validation subject, and apply the parcellations on the remaining `T-T1` sessions to compute homogeneity. We will do a grid search for `beta`, `w` and `c`, and compute homogeneity for all validation subjects. We then average the homogeneity across the validation subjects, the set of `beta`, `w` and `c` gives the highest averaged homogeneity will be selected to generate individual parcellations for the test set.

In this example, we will use the CoRR-HNU subjects 1 as the validation set. We use session 1 (i.e., sub-session 1 and 2) to estimate parcellation and use session 2 to compute the homogeneity.


In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step3_generate_ind_parcellations
```

Start Matlab, in Matlab command window, run the following commands:

```
project_dir = '<output_dir>/generate_individual_parcellations';

% Generate individual parcellation for validation subjects using 2 sub-sessions of session 1 and compute homogeneity using fMRI data of session 2
% Do a grid search for w and c for a specific beta
curr_beta = '5';
w_set = [40 50 60];
c_set = [40 50 60];
for i = 1:length(w_set)
    curr_w = num2str(w_set(i));
    for j = 1:length(c_set)
        curr_c = num2str(c_set(j));
        for sub = 1
            homo_with_weight(sub,:) = CBIG_ArealMSHBM_parameters_validation(project_dir,'fsaverage6','2','400',num2str(sub), curr_w, curr_c, curr_beta, 'gMSHBM');
        end
        homo(i,j) = mean(mean(homo_with_weight));
    end
end
 
```
The set of `w`, `c` and `beta` with the highest `homo` should be selected as the parameters for test set. If `w = 50`, `c = 50`, for subject 1, `homo_with_weight = 0.6042`.

The individual parcellation of validation set will be saved into

+ `<output_dir>/generate_individual_parcellations/ind_parcellation_gMSHBM/validation_set/2_sess/beta5/Ind_parcellation_MSHBM_sub1_w50_MRF50_beta5.mat`

The relevant homogeneity results will be saved into

+ `<output_dir>/generate_individual_parcellations/homogeneity_gMSHBM/validation_set/2_sess/beta5/Ind_homogeneity_ArealMSHBM_sub1_w50_MRF50_beta5.mat`

#### 4.2 Test set

We will still use gMSHBM as an example. Assuming `w = 50`, `c = 50` and `beta = 5` is the set of parameters give the highest homogeneity. We then apply them on the test set to generate parcellations. In this example, we use CoRR-HNU subject 2 as the test subject.

The `profile_list/test_set` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/profile_list/test_set/lh_sess?.txt`
+ `<output_dir>/generate_individual_parcellations/profile_list/test_set/rh_sess?.txt`

contain functional connectivity profiles of 2 CoRR-HNU subjects with 2 sessions. 

```
project_dir = '<output_dir>/generate_individual_parcellations';

% Generate individual parcellation for subject 2 using 4 sub-sessions
[lh_labels, rh_labels] = CBIG_ArealMSHBM_gMSHBM_generate_individual_parcellation(project_dir,'fsaverage6','4','400','2','50','50','5','test_set');

% Visualize the parcellation of subject 2
CBIG_DrawSurfaceMaps(lh_labels,rh_labels, 'fsaverage6', 'inflated', 1, 400, rand(401,3));
```

The generated individual parcellations will be saved under:
+ `<output_dir>/generate_individual_parcellations/ind_parcellation/test_set/4_sess/beta5/Ind_parcellation_MSHBM_sub2_w50_MRF50_beta5.mat`

### Check example results
----

+ We provide a script which checks user's example results with our reference results:

   `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples/CBIG_ArealMSHBM_check_example_results.m`

   Assume your example results are saved in `out_dir`. Please run the following command to check your results:
   ```
   cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/examples
   CBIG_ArealMSHBM_check_example_results(out_dir)
   ```

----

Bugs and Questions
====
Please contact Ru(by) Kong at roo.cone@gmail.com.

