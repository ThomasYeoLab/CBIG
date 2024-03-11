function [group_FC, group_FCD_CDF] = CBIG_pFIC_generate_HCP_group_level_FC_FCD(TC_path, subject_list)

% [group_FC, group_FCD_CDF] = CBIG_pFIC_generate_HCP_group_level_FC_FCD(subject_list)
% This function generates a group-average FC and group-average FCD CDF based on the given list
% of subjects
%
% Input:
%   - TC_path: absolute path to the Desikan parcellated time series (a 72-by-1200 .mat file). Refer to 
%       CBIG_generate_fslr_TC_Desikan.m for more details
%   - subject_list: absolute path to an n-by-1 vector where each row is a subject ID (n = number of subjects)
% Output:
%   - group_FC: an 68-by-68 group-level FC matrix
%   - group_FCD_CDF: a 10000-by-1 group-level FCD culumative distribution function (CDF)
%
% Example:
% [FC_train, FCD_train] = CBIG_pFIC_generate_HCP_group_level_FC_FCD('HCP_TC/Desikan', ...
%   'Zhang2023_pFIC/replication/HCP/input/training_subject_list.txt');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

FC_all = [];
FCD_CDF_all = [];
num_roi = 68;

% read subject list
fid = fopen(subject_list, 'r');
stringPattern = textscan(fid,'%s');
fclose(fid);
subject_list = string(stringPattern{:});

for i = 1:length(subject_list)
   disp(['Processing: ' num2str(i) '/' num2str(length(subject_list)) ' subjects...'])
   subject = subject_list{i}; 
   for session = {'1_LR', '1_RL', '2_LR', '2_RL'}
       input_path = [TC_path '/' subject '_rfMRI_REST' session{:} '.mat'];
       if exist(input_path, 'file')
          TC = load(input_path);
          TC = TC.TC;
          TC([1, 5, 37, 41], :) = []; % remove medial wall
          % only consider runs with 1200 frames to ensure that FCD matrix has consistent dimensions across runs
          if size(TC, 2) == 1200
              TC = TC';
              % this function can be found in our GitHub repo
              % ThomasYeoLab/CBIG/utilities/matlab/stats/CBIG_self_corr.m
              FC = CBIG_self_corr(TC); 
              FC_all = cat(3, FC_all, FC);

              FCD_matrix = zeros(num_roi*(num_roi-1)/2, 1200 - 82);
              % sliding window, slide one TR at a time
              % Each sliding window covers roughly 1-min of BOLD signals. Since HCP time series has a TR of 0.72s,
              % that is equivalent to approx. 83 frames
              for j = 1:(1200 - 82)
                    TC_section = TC(j:j+82, :);
                    FC_section = CBIG_self_corr(TC_section);
                    FC_vec_section = FC_section(triu(true(size(FC_section, 1)), 1)); 
                    FCD_matrix(:, j) = FC_vec_section;
              end

              FCD_matrix = corr(FCD_matrix); 
              % extract the upper triangle of FCD matrix (excluding the main diagonal)
              FCD_vec = FCD_matrix(triu(true(size(FCD_matrix, 1)), 1)); 

              % turn it into a PDF, 10000 bins 
              FCD_PDF = histcounts(sort(FCD_vec), -1:0.0002:1);

              % then turn it into a CDF
              FCD_CDF = cumsum(FCD_PDF); 
              FCD_CDF_all = cat(1, FCD_CDF_all, FCD_CDF);
          end
       end
   end  
end

% Do Fisher transformation before averaging FCs, then inverse Fisher
% transformation after averaing FCs
group_FC = tanh(mean(CBIG_StableAtanh(FC_all), 3));
group_FCD_CDF = mean(FCD_CDF_all);
group_FCD_CDF = round(group_FCD_CDF');

end
