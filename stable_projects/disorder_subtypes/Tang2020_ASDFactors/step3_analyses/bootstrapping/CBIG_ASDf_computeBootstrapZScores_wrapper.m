function CBIG_ASDf_computeBootstrapZScores_wrapper(in_dir, out_dir)
% CBIG_ASDf_computeBootstrapZScores_wrapper(in_dir, out_dir)
% 
% Wrapper function to compute bootstrapped z-scores and p-values
%
% Input:
%     - in_dir:
%           Absolute path to directory where bootstrapped factor estimates
%           are saved
%     - out_dir:
%           Absolute path to directory where output results will be saved
%
% Example:
%       CBIG_ASDf_computeBootstrapZScores_wrapper('~/bootstrapping/estimates', '~/bootstrapping/analyses')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

N = 100; % 100 resamples
k = '3'; % number of factors
scalelim_blk = [];
M = 419*418/2;
M_blk = 18*19/2;
is_sum = false;

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR, 'step2_polarLDA'));
addpath(fullfile(CODE_DIR, 'step3_analyses','utilities'));

%% Original factor
inputDir_best = fullfile(CODE_DIR,'data_release','files_polarLDA','model_K3');
beta_best = exp(load(fullfile(inputDir_best,'final.beta')));
rho_best = exp(load(fullfile(inputDir_best,'final.rho')));
Mean_best = beta_best.*(2*rho_best-1);
Mean_best = Mean_best([3 2 1],:); % re-order the original factors

Mean_best_net = zeros(3,M_blk);
Mean_best_net(1,:) = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly(Mean_best(1,:), scalelim_blk, [], is_sum);
Mean_best_net(2,:) = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly(Mean_best(2,:), scalelim_blk, [], is_sum);
Mean_best_net(3,:) = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly(Mean_best(3,:), scalelim_blk, [], is_sum);
close all;

f1_all = zeros(N, M);
f2_all = zeros(N, M);
f3_all = zeros(N, M);

f1_net = zeros(N, M_blk);
f2_net = zeros(N, M_blk);
f3_net = zeros(N, M_blk);

% Loop over all bootstrapped samples
for i = 1:N
    input_dir = fullfile(in_dir, ['resampled_' num2str(i)]);
    
    r_dir = 'init_with_model';
    
    beta = exp(load(fullfile(input_dir,sprintf('k%s/%s', k, r_dir),'final.beta')));
    rho = exp(load(fullfile(input_dir,sprintf('k%s/%s', k, r_dir),'final.rho')));
    Mean = beta.*(2*rho-1);
    
    order = CBIG_ASDf_hunMatch(str2double(k), Mean_best, Mean);
    Mean = Mean(order,:);
    
    % save current bootstrapped sample factors
    f1_all(i,:) = Mean(1,:);
    f2_all(i,:) = Mean(2,:);
    f3_all(i,:) = Mean(3,:);
    
    % 18x18 blocks
    f1_net(i,:) = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly(Mean(1,:), scalelim_blk, [], is_sum);
    f2_net(i,:) = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly(Mean(2,:), scalelim_blk, [], is_sum);
    f3_net(i,:) = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_NetworksOnly(Mean(3,:), scalelim_blk, [], is_sum);
    close all;
end

save(fullfile(out_dir,'bootstrappedFactors.mat'), 'f1_all', 'f2_all', 'f3_all');
save(fullfile(out_dir,'bootstrappedFactorsByNetworks.mat'), 'f1_net', 'f2_net', 'f3_net');

%% 419x419 mean and std
means = zeros(3,M);
stds = zeros(3,M);
means(1,:) = mean(f1_all,1);
means(2,:) = mean(f2_all,1);
means(3,:) = mean(f3_all,1);

stds(1,:) = std(f1_all, 0, 1);
stds(2,:) = std(f2_all, 0, 1);
stds(3,:) = std(f3_all, 0, 1);

%% 18x18 mean and std
means_blk = zeros(3, M_blk);
stds_blk = zeros(3, M_blk);
means_blk(1,:) = mean(f1_net,1);
means_blk(2,:) = mean(f2_net,1);
means_blk(3,:) = mean(f3_net,1);

stds_blk(1,:) = std(f1_net, 0, 1);
stds_blk(2,:) = std(f2_net, 0, 1);
stds_blk(3,:) = std(f3_net, 0, 1);

save(fullfile(out_dir,'meanStd.mat'), 'means', 'stds');
save(fullfile(out_dir,'meanStdByNetworks.mat'), 'means_blk', 'stds_blk');

%% Compute z-scores
z_scores = Mean_best ./ stds;
z_scores_blk = Mean_best_net ./ stds_blk;

% p-value, two-tailed
p_vals = 2 * (1 - normcdf(abs(z_scores)));
save(fullfile(out_dir,'bootstrapped_pVals.mat'), 'z_scores', 'p_vals');

p_vals_blk = 2 * (1 - normcdf(abs(z_scores_blk)));
save(fullfile(out_dir,'bootstrappedByNetworks_pVals.mat'), 'z_scores_blk', 'p_vals_blk');

close all;

%% Remove paths
rmpath(fullfile(CODE_DIR, 'step2_polarLDA'));
rmpath(fullfile(CODE_DIR, 'step3_analyses','utilities'));
