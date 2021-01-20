function CBIG_MMLDA_example_wrapper(out_dir, queue)
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])


doc_dir = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/examples/input'];
ref_dir = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/examples/correct_output'];
%%%
% 1. MMLDA estimation
%%%
mmlda_dir = [out_dir '/estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];
mkdir(mmlda_dir)
cmd=['${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/' ...
    'step2_MMLDA/CBIG_MMLDA_est.sh ' ...
    '-a ' doc_dir '/ADNI2_bl_AD_meanCNstdALL_plus1_brain_2sub.dat ' ...
    '-b ' doc_dir '/ADNI2_bl_AD_meanCNstdALL_plus1_behavior_2sub.dat ' ...
    '-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt ' ...
    '-k 2 ' ...
    '-s 1 ' ...
    '-e 5 ' ...
    '-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA ' ...
    '-o ' out_dir '/estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub '];
if exist('queue', 'var')
    cmd = [cmd '-q ' queue];
end
system(cmd)
cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n MMLDA_est';
system(cmd)

%%%
% 2. Visualize factors
%%%
t = readtable([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
'Sun2019_ADJointFactors/utilities/Behavior_Abbreviation.csv']);
behavior_name = t.Behavioral_Name_Short;
behavior_domain = t.Behavioral_Cate;

% ADNI2 bl K=2 2 subjects
visualize_dir = [out_dir '/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];
mkdir(visualize_dir)
mask = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/' ...
'ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz'];
min_thresh = 7.5e-6;
max_thresh = 1.5e-5;
k = 2;
CBIG_MMLDA_visualize_factors(mmlda_dir, visualize_dir, k, min_thresh, max_thresh, ...
 mask, behavior_name, behavior_domain)

%%%
% 3. MMLDA inference
%%% 
inf_dir = [out_dir '/inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];
rname = ls([visualize_dir '/k2/']);
cmd=['${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/' ...
    'step2_MMLDA/CBIG_MMLDA_inf.sh ' ...
    '-a ' doc_dir '/ADNI2_bl_AD_meanCNstdALL_plus1_brain_1sub.dat ' ...
    '-b ' doc_dir '/ADNI2_bl_AD_meanCNstdALL_plus1_behavior_1sub.dat ' ...
    '-t ${CBIG_CODE_DIR}/external_packages/mmlda-c-dist/code/settings-100iter.txt ' ...
    '-k 2 ' ...
    '-d ' visualize_dir '/k2/' rname(1:2) '/final ' ...
    '-m ${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA ' ...
    '-o ' out_dir '/inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub ' ...
    '-n ADNI2_bl_AD_meanCNstdALL_plus1_1sub '];
if exist('queue', 'var')
    cmd = [cmd '-q ' queue];
end
system(cmd)
cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n MMLDA_inf';
system(cmd)
%%%
% 4. Compare output results with reference
%%%
% Check if final estimated model parameters are the same as reference
fprintf('\n--------[CHECK 1] Checking if estimated model parameters are the same as reference:\n');

curr_out_dir = [out_dir '/estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3'];
curr_ref_dir = [ref_dir '/estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3'];

% beta1 estimate
if (~exist([curr_out_dir '/final.beta1'], 'file'))
    fprintf('[FAILED] Final beta1 estimate result missing!\n');
else
    beta1 = load([curr_out_dir '/final.beta1']);
    ref_beta1 = load([curr_ref_dir '/final.beta1']);
    rel_diff_beta1 = max(max(abs((beta1 - ref_beta1) ./ ref_beta1)));
    if (rel_diff_beta1 > 1e-4)
        fprintf(['[FAILED] Final beta1 estimate is different from reference,' ...
            ' max abs relative diff: %f.\n'], rel_diff_beta1);
    else
        fprintf('[PASSED] Final beta1 estimate is the same as reference.\n');
    end
    % release some memory
    clear beta1;
    clear ref_beta1;
end

% beta2 estimate
if (~exist([curr_out_dir '/final.beta2'], 'file'))
    fprintf('[FAILED] Final beta2 estimate result missing!\n');
else
    beta2 = load([curr_out_dir '/final.beta2']);
    ref_beta2 = load([curr_ref_dir '/final.beta2']);
    rel_diff_beta2 = max(max(abs((beta2 - ref_beta2) ./ ref_beta2)));
    if (rel_diff_beta2 > 1e-4)
        fprintf(['[FAILED] Final beta2 estimate is different from reference,' ...
            ' max abs diff: %f.\n'], rel_diff_beta2);
    else
        fprintf('[PASSED] Final beta2 estimate is the same as reference.\n');
    end
    % release some memory
    clear beta2;
    clear ref_beta2;
end

% gamma estimate
if (~exist([curr_out_dir '/final.gamma'], 'file'))
    fprintf('[FAILED] Final gamma estimate result missing!\n');
else
    gamma = load([curr_out_dir '/final.gamma']);
    ref_gamma = load([curr_ref_dir '/final.gamma']);
    rel_diff_gamma = max(max(abs((gamma - ref_gamma) ./ ref_gamma)));
    if (rel_diff_gamma > 1e-4)
        fprintf(['[FAILED] Final gamma estimate is different from reference,' ...
            ' max abs diff: %f.\n'], rel_diff_gamma);
    else
        fprintf('[PASSED] Final gamma estimate is the same as reference.\n');
    end
    % release some memory
    clear gamma;
    clear ref_gamma;
end
        
% likelihood
if (~exist([curr_out_dir '/likelihood.dat'], 'file'))
    fprintf('[FAILED] Likelihood file missing!\n');
else
    likelihood = load([curr_out_dir '/likelihood.dat']);
    likelihood = likelihood(:, 1);
    ref_likelihood = load([curr_ref_dir '/likelihood.dat']);
    ref_likelihood = ref_likelihood(:, 1);
    rel_diff_likelihood = max(max(abs((likelihood - ref_likelihood) ./ ref_likelihood)));
    if (rel_diff_likelihood > 1e-4)
        fprintf(['[FAILED] Likelihood is different from reference,' ...
            ' max abs diff: %f.\n'], rel_diff_likelihood);
    else
        fprintf('[PASSED] Likelihood is the same as reference.\n');
    end
    % release some memory
    clear likelihood;
    clear ref_likelihood;
end

% Check if Pr(Factor|Participant) are the same as reference
fprintf('\n--------[CHECK 2] Checking if Pr(Factor|Participant) are the same:\n');

curr_out_dir = [out_dir '/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3'];
curr_ref_dir = [ref_dir '/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3'];

factorComp = load([curr_out_dir '/FactorComp.txt']);
ref_factorComp = load([curr_ref_dir '/FactorComp.txt']);
rel_diff_factorComp = max(max(abs((factorComp - ref_factorComp) ./ ref_factorComp)));
if (rel_diff_factorComp > 1e-4)
    fprintf(['[FAILED] Factor composition is different from reference,' ...
        ' max abs diff = %f.\n'], rel_diff_factorComp);
else
    fprintf('[PASSED] Factor composition is the same as reference.\n');
end
% release some memory
clear factorComp;
clear ref_factorComp;

% Check if inferred gamma file is the same as reference
fprintf('\n--------[CHECK 3] Checking if inferred factor composition is the same:\n');

curr_out_dir = [out_dir '/inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];
curr_ref_dir = [ref_dir '/inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];

% Gamma estimated by joint modalities
if (~exist([curr_out_dir '/k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat'], 'file'))
    fprintf('[FAILED] Final gamma inference result missing!\n');
else
    gamma = load([curr_out_dir '/k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat']);
    ref_gamma = load([curr_ref_dir '/k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat']);
    rel_diff_gamma = max(max(abs((gamma - ref_gamma) ./ ref_gamma)));
    if (rel_diff_gamma > 1e-4)
        fprintf(['[FAILED] Final gamma inference is different from reference,' ...
            ' max abs diff: %f.\n'], rel_diff_gamma);
    else
        fprintf('[PASSED] Final gamma inference is the same as reference.\n');
    end
    % release some memory
    clear gamma;
    clear ref_gamma;
end

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])


