function CBIG_TRBPC_compute_PFM_average(PFM_dir, outdir, hypothesis_ind, datadriven_ind)
% CBIG_TRBPC_compute_PFM_average(PFM_dir, outdir, hypothesis_ind, datadriven_ind)
% 
% This function prepares the files that will be used for the feature transfer
% learning in Chen & Tam 2021 paper. For each behavior, we compute the
% average PFM of all behaviors in a behavior domain excluding the behavior itself.
%
% Inputs:
%   - PFM_dir
%     Directory where regression and PFM results are stored
%
%   - outdir
%     output directory for the PFM average files
%
%   - hypothesis_ind
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for hypothesis-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
%   - datadriven_ind
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for data-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
% Outputs:
%   Output files will be saved in outdir
%
% Written by Jianzhong Chen & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

N_edge = 419*418/2;
PFM_all_score = zeros(N_edge*4,120,36);
for i = 1:36
    load(fullfile(PFM_dir,['PFM_score' num2str(i) '_all_folds.mat']), 'PFM_all_folds');
    PFM_all_score(:,:,i) = PFM_all_folds;
end

%% hypothesis-driven
PFM_all_folds = zeros(N_edge*4,120,36);
cat = {'cog','mh','pers'};
for c = 1:3
    disp(cat{c});
    for i = 1:36
        behav_ind = setdiff(hypothesis_ind.(cat{c}), i);
        mean_rel = mean(PFM_all_score(:,:,behav_ind),3);
        PFM_all_folds(:,:,i) = squeeze(mean_rel);
    end
    save(fullfile(outdir, ['PFM_domain_average_hypothesis_all_folds_' cat{c} '.mat']), 'PFM_all_folds', '-v7.3');
end

%% data-driven
PFM_all_folds = zeros(N_edge*4,120,36);
cat = {'cog','mh','pers'};
for c = 1:3
    disp(cat{c});
    for i = 1:36
        behav_ind = setdiff(datadriven_ind.(cat{c}), i);
        mean_rel = mean(PFM_all_score(:,:,behav_ind),3);
        PFM_all_folds(:,:,i) = squeeze(mean_rel);
    end
    save(fullfile(outdir, ['PFM_domain_average_datadriven_all_folds_' cat{c} '.mat']), 'PFM_all_folds', '-v7.3');
end
