function CBIG_ICCW_LRR_get_parameters(KRR_dir, outdir, behav_ind)

% function CBIG_ICCW_LRR_get_parameters(KRR_dir, outdir, behav_ind)
%
% This function prepares the input parameters for LRR. The script
% copies the param file from KRR and adds the additional arguments
% required to run LRR.
%
% Note that LRR is run using Elasticnet with alpha set to 0,
% therefore only L2 regularization is used.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% load existing param file
load(fullfile(KRR_dir, 'param.mat'),'param');
% load FC matrix
load(fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', ...
    'ChenTam2022_TRBPC', 'FC', 'FC_subjects_rs_all_score_mf.mat'));
% normalize FC matrix
FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);
% add to param file
param.feature_mat = FC_all;
param.y = param.y(:,behav_ind);
param.outdir=fullfile(outdir,['behav_' num2str(behav_ind)]);
param.split_name = 'rng1';
param.alpha = 0;
mkdirp(param.outdir)
% save new param file
save(fullfile(param.outdir,'param.mat'),'param','-v7.3');

end