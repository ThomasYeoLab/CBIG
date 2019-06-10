# Infer factor compositions of new subjects
Even though we provide the factor compositions of ADNI2 bl/m12 and ADNI1 bl subjects, the user
may want to infer factor compositions of an unseen subject by using our estimated factors. 

----

# Usage
The main function `CBIG_MMLDA_infer_new.sh` will do the following things
1. SPM VBM preprocessing of new subjects by using our ADNI2 bl custom template.
2. Convert grey matter probability maps and behavioral scores to documents by 
   using our ADNI2 bl CN parameters.
3. Infer factor compositions of new subjects by using our factors.
4. Normalize the gamma file into probabilites and reorder it to match our paper.

The user can check `CBIG_MMLDA_infer_new_subjects_wrapper.sh` for an example usage.
If the user wants to estimate new factors and do inference on new factors,
please check `step2_MMLDA` folder.

### Input
There are three input files:
1. `T1_list.txt`:   Text file with each line being the path to a T1 image
2. `id_list.txt`:   Text file with each line being the id of the subject
3. `subinfo.csv`:   csv file which includes subject information, such as id, age, sex, behavior scores
                    For the id column, it mush match the idList text file.
                    For the format of the behavior scores, please check the `input/subinfo.csv`
                    You have to fill in all behavior scores except the 'ADAS: Recall instruction, ANART, CFT: vegetable'
                    columns. Because we do not use these three columns in following processing, you can just fill in these three columns with NaN.

For the format of input files, the user can check `input` folder for an example.

In addition, the user can run `./CBIG_MMLDA_infer_new_subjects.sh` to see more details.

### Output
In the output folder, the `prob` folder will contain following files
1. `k3_inf_mmlda_prob.csv`: joint factor loadings by using both atrophy and behavior
2. `k3_inf1_mmlda_prob.csv`: factor loading by using atrophy
3. `k3_inf2_mmlda_prob.csv`: factor loading by using behavior
In the csv files, there are three columns corresponding to MTL-Memory, 
Lateral Temporal-Language, Posterior Cortical-Executive factor compositions respectively.

### Missing data
If the user have few missing data in the behavioral scores, the user can impute the scores by using
GLM or other methods. The matlab function for imputing behavioral scores using GLM can be found in 
`Sun2019_ADJointFactors/utilities/CBIG_MMLDA_matrix_completion_GLM.m`.  
If the user have a lot of missing data in the behavioral scores, the user
can infer factor loadings by only using atrophy maps. This means the user can fill in arbitrary
values in the behavioral scores and only check the `k3_inf1_mmlda_prob.csv` for factor loadings.
