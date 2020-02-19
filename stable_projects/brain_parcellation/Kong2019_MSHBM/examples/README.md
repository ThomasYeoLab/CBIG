# Examples of Multi-session hierarchical Bayesian model (MS-HBM)

The examples of MS-HBM contain three tasks:
1. Generating profiles and estimating initialization parameters
   - CoRR-HNU data in fsaverage5 surface space: 2 subjects, each subject has 2 sessions, each subject has 1 run
2. The group priors estimation on example data:
   - CoRR-HNU data in fsaverage5 surface space: 2 subjects, each subject has 2 sessions
3. The individual-level parcellation generation on example data:
   - CoRR-HNU data in fsaverage5 surface space: 2 subjects, each subject has 2 sessions
   
Note that for simplicity, the example data of task 1 and task 2 are the same. However, the initialization parameters estimated in task 1 **won't** be used in task 2, the group priors estimated in task 2 **won't** be used in task 3. We will use the initialization parameters and group priors pre-computed by using GSP dataset (37 subjects, 2 sessions). 

----

References
==========
+ Kong R, Li J, Orban C, et al. [Spatial Topography of Individual-Specific Cortical Networks Predicts Human Cognition, Personality, and Emotion](https://academic.oup.com/cercor/article/29/6/2533/5033556). Cerebral Cortex, 29(6):2533-2551, 2019

----

Data
====
In this example, we will perform the algorithms on the functional connectivity profiles of CoRR-HNU datasets:
+ **CoRR-HNU:**
  
  `$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0?/subj0?_sess?`

----

Run
====

### Generating input data
----

In the terminal, specify the output directory and call the example data generation script:
```
$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/examples/CBIG_MSHBM_create_example_input_data.sh <output_dir>
```
Note that in your `<output_dir>`, there will be three folders:
+ `<output_dir>/generate_profiles_and_ini_params`
+ `<output_dir>/estimate_group_priors`
+ `<output_dir>/generate_individual_parcellations`

We will then provide examples for each task separately.

### Generating profiles and initialization parameters
----

After [**Generating input data**](#generating-input-data),

the `data_list` folder under `<output_dir>/generate_profiles_and_ini_params`:

+ `<output_dir>/generate_profiles_and_ini_params/data_list/fMRI_list/?h_sub?_sess?.txt`

contains the left and right hemisphere fMRI data in `fsaverage5` for each subject and each session. Note that our example data only has 1 run for each session, if there are multiple runs, then the corresponding `?h_sub?_sess?.txt` should have multiple columns, where each column corresponds to each run, the columns should be separated by white space.

+ `<output_dir>/generate_profiles_and_ini_params/data_list/censor_list/sub?_sess?.txt`

contains the censor lists for each subject and each session. The outlier frames will be denoted as 0s. Again, if there are multiple runs, `sub?_sess?.txt` should have multiple columns. Furthermore, if the user don't want to mask out the outliers in computing profiles, just leave the `censor_list` folder empty.

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/step1_generate_profiles_and_ini_params
```

Start Matlab, in Matlab command window, run the following commands to generate profiles for each session of the 2 subjects, the fMRI data are in `fsaverage5`, the ROIs are defined by `fsaverage3`.

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
for sub = 1:2
 for sess = 1:2
	CBIG_MSHBM_generate_profiles('fsaverage3','fsaverage5',project_dir,num2str(sub),num2str(sess),'0');
 end
end
```
The profiles will be saved into `profiles` folder:

+ `<output_dir>/generate_profiles_and_ini_params/profiles/sub?/sess?/?h.sub?_sess?_fsaverage5_roifsaverage3.surf2surf_profile.nii.gz`


In practice, if the user only has 1 session and would like to apply it on MSHBM, the user need to split the session into two sub-sessions to make MSHBM work. In this example, the CoRR-HNU subjects will also be utilized as the example data of validation set, we will use session 1 of subjects 1 and 2 to generate parcellations. Note that in practice, the validation and test set should be **different** sets.   

To do that, run the following commands to split the first session of the 2 CoRR-HNU subjects and generate 2 profiles for each sub-session.

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
for sub = 1:2
 CBIG_MSHBM_generate_profiles('fsaverage3','fsaverage5',project_dir,num2str(sub),'1','1');
end
```
The profiles will be saved into `profiles` folder:

+ `<output_dir>/generate_profiles_and_ini_params/profiles/sub?/sess?/?h.sub?_sess?_fsaverage5_roifsaverage3.surf2surf_profile_1.nii.gz`

+ `<output_dir>/generate_profiles_and_ini_params/profiles/sub?/sess?/?h.sub?_sess?_fsaverage5_roifsaverage3.surf2surf_profile_2.nii.gz`

To generate the initialization parameters, we will apply Yeo2011 clustering algorithm on group averaged profiles. Run the following command to generate averaged profiles across 2 CoRR-HNU subjects with 2 sessions. Note that the profiles of split sub-sessions will be ignored. In practice, if the subject has different number of sessions, please set the `num_sess` to be the maximum number of sessions.

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
num_sub = '2';
num_sess = '2';
CBIG_MSHBM_avg_profiles('fsaverage3','fsaverage5',project_dir,num_sub,num_sess);
```

The results will be saved into `profiles/avg_profile` folder:

+ `<output_dir>/generate_profiles_and_ini_params/profiles/avg_profile/?h_fsaverage5_roifsaverage3_avg_profile.nii.gz`

We will then apply Yeo2011 clustering algorithm on above averaged profiles to generate the initialization parameters with 17 networks and 2 random initialization:

```
project_dir = '<output_dir>/generate_profiles_and_ini_params';
num_clusters = '17';
num_initialization = '2';
CBIG_MSHBM_generate_ini_params('fsaverage3','fsaverage5',num_clusters,num_initialization, project_dir)
```

The results will be saved into `group` folder:

+ `<output_dir>/generate_profiles_and_ini_params/group/group.mat`

contains the estimated initialization parameter `clustered.mtc`, which is the group-level functional connectivity profile of networks. `group.mat` also contains the estimated group-level parcellation `lh_labels` and `rh_labels`. The visualization of group-level parcellation is:

![visualization_of_group](results/generate_profiles_and_ini_params/group/group_figure.jpg)

Note that there are only 2 random initializations and averaged profiles of 2 subjects in this example, so the clustering results are not optimal. In practice, we suggest the user to at least run 1000 random initializations. Please also note that we **will not** use the results generated in this step in the following sections.

### Group priors estimation
----

Note that the functional connectivity profiles and the initialization parameters from the Yeo2011 group clustering are already generated. There is no need to re-generate the profiles. 

After [**Generating input data**](#generating-input-data),

the `group` folder under `<output_dir>/estimate_group_priors`:
+ `<output_dir>/estimate_group_priors/group/group.mat`

contains the **pre-computed** initialization parameters from group-level 17-network clustering generated by Yeo et al., 2011.

The `profile_list/training_set` folder under `<output_dir>/estimate_group_priors`:
+ `<output_dir>/estimate_group_priors/profile_list/training_set/lh_sess?.txt`
+ `<output_dir>/estimate_group_priors/profile_list/training_set/rh_sess?.txt`

contain functional connectivity profiles of 2 CoRR-HNU subjects with 2 sessions. These files are also **pre-generated**.

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/step2_estimate_priors
```

Start Matlab, in Matlab command window, run the following commands to estimate group priors on 2 subjects with 2 sessions for 17 networks, we will set the maximum number of iteration to be 5:

```
project_dir = '<output_dir>/estimate_group_priors';
Params = CBIG_MSHBM_estimate_group_priors(project_dir,'fsaverage5','2','2','17','5');
```
As there are only two subjects, the script may not be able to converge, we will only run 5 iterations in this example. The results of each iteration will be saved into `priors` folder, the results of each iteration will be saved as `Params_iteration?.mat`, which contains a struct variable `Params`. The final estimated group priors should be saved as `Params_Final.mat`.

+ `<output_dir>/estimate_group_priors/priors`

The estimated group priors include:
1) Inter-subject functional connectivity variability -- `Params.epsil`
2) Group-level connectivity profiles for each network -- `Params.mu`
3) Intra-subject functional connectivity variability -- `Params.sigma`
4) Spatial prior which denotes the probability of each network occurring at each location -- `Params.theta`

The user can compare their estimated `Params_iteration?.mat` with 
+ `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/examples/results/estimate_group_priors/priors/Params_iteration?.mat`

### Individual-level parcellations generation
----

Note that the group priors estimated in **Group priors estimation** won't be used to generate individual parcellations. We will use the group priors pre-computed by using GSP dataset (37 subjects, 2 sessions), which are saved in:

+ `<output_dir>/generate_individual_parcellations/priors/Params_Final.mat`

#### Validation set

This algorithm contains two parameters:
+ `w`: The weight of group spatial prior `Params.theta`.
+ `c`: The weight of MRF smoothness prior.

These two parameters can be selected by validation set. Assuming subjects in validation set has `T` sessions, we normally use `T1` sessions to generate the individual parcellations for each validation subject, and apply the parcellations on the remaining `T-T1` sessions to compute homogeneity. We will do a grid search for `w` and `c`, and compute homogeneity for all validation subjects. We then average the homogeneity across the validation subjects, the pair of `w` and `c` gives the highest averaged homogeneity will be selected to generate individual parcellations for the test set.

In this example, we will use the CoRR-HNU subjects 1 and 2 as the validation set. Note that even though in this example, the validation set and test set is the same, in practice, they should be **different** sets of subjects. Here, we will use session 1 of subject 1 and 2 to generate parcellations, and apply them on session 2 to compute the homogeneity.

After [**Generating input data**](#generating-input-data),

the `data_list` folder under `<output_dir>/generate_individual_parcellations`:

+ `<output_dir>/generate_individual_parcellations/data_list/validation_fMRI_list/?h_sub?.txt`

contains the left and right hemisphere fMRI data of session 2 in `fsaverage5` for each subject. Note that in this example, our validation fMRI data only has 1 session with 1 run, if there are multiple runs and multiple sessions, then the corresponding `?h_sub?.txt` should have multiple columns and rows, where each row corresponds to each session, each column corresponds to each run, the columns should be separated by white space. For example:
```
subject1_rs_fMRI_sess1_run1.nii.gz subject1_rs_fMRI_sess1_run2.nii.gz
subject1_rs_fMRI_sess2_run1.nii.gz
subject1_rs_fMRI_sess3_run1.nii.gz subject1_rs_fMRI_sess3_run2.nii.gz
```

The `profile_list/validation_set` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/profile_list/validation_set/lh_sess?.txt`
+ `<output_dir>/generate_individual_parcellations/profile_list/validation_set/rh_sess?.txt`

contain functional connectivity profiles of 2 CoRR-HNU subjects with 2 sub-sessions which are generated by splitting session 1 into two. 

In the terminal:
```
cd $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/step3_generate_ind_parcellations
```

Start Matlab, in Matlab command window, run the following commands:

```
project_dir = '<output_dir>/generate_individual_parcellations';

% Generate individual parcellation for validation subjects using 2 sub-sessions of session 1 and compute homogeneity using fMRI data of session 2
% Do a grid search for w and c
w_set = [60 80 100 120];
c_set = [30 40 50 60];
for i = 1:length(w_set)
    for j = 1:length(c_set)
        for sub = 1:2
		    homo_with_weight(sub,:) = CBIG_MSHBM_parameters_validation(project_dir,'fsaverage5','2','17',num2str(sub), num2str(w_set(i)),num2str(c_set(j)));
        end
        homo(i,j) = mean(mean(homo_with_weight));
    end
end
 
```
The pair of `w` and `c` with the highest `homo` will be selected as the parameters for test set. If `w = 100`, `c = 50`, for subject 1, `homo_with_weight = 0.3316`.


#### Test set

Assuming `w = 100` and `c = 50` is the pair of parameters give the highest homogeneity. We then apply them on the test set to generate parcellations.

The `profile_list/test_set` folder under `<output_dir>/generate_individual_parcellations`:
+ `<output_dir>/generate_individual_parcellations/profile_list/test_set/lh_sess?.txt`
+ `<output_dir>/generate_individual_parcellations/profile_list/test_set/rh_sess?.txt`

contain functional connectivity profiles of 2 CoRR-HNU subjects with 2 sessions. 

```
project_dir = '<output_dir>/generate_individual_parcellations';

% Generate individual parcellation for subject 1 using 2 sessions
[lh_labels1, rh_labels1] = CBIG_MSHBM_generate_individual_parcellation(project_dir,'fsaverage5','2','17','1','100','50');

% Generate individual parcellation for subject 2 using 2 sessions
[lh_labels2, rh_labels2] = CBIG_MSHBM_generate_individual_parcellation(project_dir,'fsaverage5','2','17','2','100','50');

% Visualize the parcellation of subject 1 and 2
group = load(fullfile(getenv('CBIG_CODE_DIR'),'/stable_projects/brain_parcellation/Kong2019_MSHBM/examples/results/estimate_group_priors/group/group.mat'));
CBIG_DrawSurfaceMaps(lh_labels1,rh_labels1, 'fsaverage5', 'inflated',0,17,group.colors);
CBIG_DrawSurfaceMaps(lh_labels2,rh_labels2, 'fsaverage5', 'inflated',0,17,group.colors);
```
The generated individual parcellations will be saved under:
+ `<output_dir>/generate_individual_parcellations/ind_parcellation`

The user can compare the generated individual parcellations with the following figures:

##### CoRR-HNU subject 1
![visualization_of_subject_1](results/generate_individual_parcellations/figures/ind_parcellation_sub1.jpg)

##### CoRR-HNU subject 2
![visualization_of_subject_2](results/generate_individual_parcellations/figures/ind_parcellation_sub2.jpg)

----

Bugs and Questions
====
Please contact Ru(by) Kong at roo.cone@gmail.com.

