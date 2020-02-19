function CBIG_ASDf_getBootstrappedSamples(output_dir, N)
% CBIG_ASDf_getBootstrappedSamples(output_dir, N)
% 
% Wrapper function to perform bootstrapping on ABIDE-II+GENDAAR ASD
% participants.
%
% Input:
%     - output_dir:
%            Absolute path to output directory where the results will be
%            saved
%     - N:
%            Interger, number of bootstrapped resamples
% 
% Example:
%       CBIG_ASDf_getBootstrappedSamples_wrapper('~/Temporary/bootstrapping', 100)
% 
%  Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Define path
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
sub_info = fullfile(UNIT_TEST_DIR,'data','subInfo_654.csv');
zScore_dir = fullfile(UNIT_TEST_DIR,'results','FC2doc','step1_output_zScores.mat');

% Add path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

% Set random seed
rng('default');

% Get diagnosis info
[~, id_dx] = CBIG_ASDf_getSubData(sub_info);
dx = cell2mat(id_dx(:,2));

% Load z scores
load(zScore_dir); % variables: z & discretized_z
z_discret = discretized_z(dx==1,:);
z = z(dx==1,:);
n_sub = size(z_discret,1);

for i = 1:N
    fprintf('====%d-th Resamples====\n', i);
    indices = randsample(n_sub, n_sub, true); % re-sample with replacement
    z_disc_resampled = z_discret(indices,:);
    z_resampled = z(indices,:);
    curr_dir = [output_dir '/resampled_' num2str(i)];
    mkdir(curr_dir);
    save([curr_dir '/zScores.mat'], 'z_disc_resampled', 'z_resampled');
    
    % Write into doc
    fileID = fopen([curr_dir '/dx1.dat'], 'w'); % clear contents
    fileID = fopen([curr_dir '/dx1.dat'], 'a'); % start appending
    for idx1 = 1:n_sub         
        fprintf('Subject order: %d \n',idx1);
        tc_one_sub = z_disc_resampled(idx1, :);
        no_terms = sum(tc_one_sub~=0);
        fprintf(fileID, '%i ', no_terms);
        for idx2 = 1:numel(tc_one_sub)
            if tc_one_sub(idx2) ~= 0
                fprintf(fileID, '%i:%i ', idx2-1, tc_one_sub(idx2));
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end

% Remove path
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
