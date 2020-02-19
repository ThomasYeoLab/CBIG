function p_tau = CBIG_MMLDA_association_factor_tau_wrapper(out_dir)
% CBIG_MMLDA_association_factor_tau_wrapper(out_dir)
%
% Investigate association between factor loadings and tau loadings.
%
% Input:
%   - out_dir   : output directory. The output will be some forest plots.
%
% Example:
%   CBIG_MMLDA_association_factor_tau_wrapper('~/example')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(genpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities']))

proj_dir = [getenv('CBIG_TESTDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];
CBIG_Sun2019_PET_preprocess = [CBIG_MMLDA_ADNI_DIR '/Sun2019_PET_preprocess/preprocessing'];

%%%%
% Get rid by considering subjects with ADNI2 preprocessed + ADNI3 original T1 
% and exclude subjects with ADNI2 original T1
%%%%
rid_257 = load([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors' ...
    '/step4_analyses/association_factor_tau/RID_ADNI23_257.csv']);

%%%%
% Get the PET tau
%%%%
% load PET tau rid list
rid_vis = importdata([CBIG_Sun2019_PET_preprocess '/lists/RID_Viscode_269all_rm_recon_fail.txt']);
rid_PET = strtok(rid_vis, '_');
rid_PET = cellfun(@str2num, rid_PET);

% load factors
beta1 = exp(load([proj_dir '/step2_MMLDA/results/results_long/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1/k3/r9/final.beta1']));
beta1_order = beta1([1 3 2], :);

% load mask
mask_file = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/ADNIDataRelease' ...
    '/SPM_VBM_files/ADNI2_bl_gm_mask_MNI2mm.nii.gz'];
mri = MRIread(mask_file);
mask = logical(mri.vol(:));

% get PET tau
mri = MRIread([CBIG_Sun2019_PET_preprocess '/output/merg_SPMVBMdeform.nii.gz']);
vol = mri.vol;
volsize = size(vol);
PETtau2d = reshape(vol, prod(volsize(1:3)), volsize(4));
PETtau2d_mask = PETtau2d(mask, :);
total_tau = sum(PETtau2d_mask)';

% Get top 5% PET tau
percent = 95;
for t = 1:3
    thr_t = prctile(beta1_order(t, :), percent);
	ind_t = beta1_order(t, :) > thr_t;
	mean_tau_prc(:, t) = sum(PETtau2d_mask(ind_t, :))' / sum(ind_t); 
end

%%%%
% Get atrophy and behavioral factor loadings
%%%%
proj_dir = [getenv('CBIG_TESTDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];
order = [1 3 2];

rid_file = [proj_dir '/step2_MMLDA/results/results_long/BrainBehavior2doc/ADNI2_PETtau_ALL_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/results_long/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf_ADNI2_PETtau_ALL_meanCNstdALL_plus1-gamma.dat'];
rid_prob_PETtau_inf = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

%%%%
% Get the intersection between rid PET and rid factor
%%%%
rid_int = intersect(rid_257, intersect(rid_PET, rid_prob_PETtau_inf(:, 1)));
ind_PET = CBIG_MMLDA_find_array_in_array(rid_int, rid_PET);
total_tau_int = total_tau(ind_PET, :);
mean_tau_prc_int = mean_tau_prc(ind_PET, :);
ind_prob = CBIG_MMLDA_find_array_in_array(rid_int, rid_prob_PETtau_inf(:, 1));
rid_prob_PETtau_inf_int = rid_prob_PETtau_inf(ind_prob, :);

%%%
% Get age, sex, dx from subinfo
%%%
subinfo_dir = [getenv('CBIG_TESTDATA_DIR') '/stable_projects/disorder_subtypes' ...
'/Sun2019_ADJointFactors/step2_MMLDA/data'];
subinfo = csvread([subinfo_dir '/ADNI23_PETtau_subinfo.csv'], 1, 0);
ind = CBIG_MMLDA_find_array_in_array(rid_int, subinfo(:, 1));
age = subinfo(ind, 3);
sex = subinfo(ind, 2);
dx = subinfo(ind, 4);

%%%%
% GLM between PET tau and atrophy/behavioral factor loadings
%%%%
k = 3;
X = [mean_tau_prc_int(dx==2, :) age(dx==2, :) sex(dx==2, :) total_tau_int(dx==2, :)];
dependency = 'independent';
compare_factor_names = {'L-M', 'C-M', 'C-L'};
CI_scale = 1;
p_thresh = 0.0167;

y = rid_prob_PETtau_inf_int(dx==2, 2);
out_name = 'association_MTL-Memory_tau';
p_mem = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, out_name, compare_factor_names, CI_scale, p_thresh)

y = rid_prob_PETtau_inf_int(dx==2, 3);
out_name = 'association_LT-Language_tau';
p_lan = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, out_name, compare_factor_names, CI_scale, p_thresh)

y = rid_prob_PETtau_inf_int(dx==2, 4);
out_name = 'association_PC-Executive_tau';
p_ef = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, out_name, compare_factor_names, CI_scale, p_thresh)

p_tau = [p_mem; p_lan; p_ef];

rmpath(genpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities']))
