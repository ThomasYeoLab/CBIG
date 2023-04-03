function CBIG_hMRF_check_replication_results(output_dir)
% CBIG_hMRF_check_replication_results(output_dir)
%
% This is the wrapper function to set up and run the clustering algorithm based on the input arguments.
% 
% Input
%   - output_dir: (string) 
%     ABSOLUTE path to the directory to which the output results will be saved.
%
% Example
%   - CBIG_hMRF_check_replication_results(your_replication_output_dir)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
if(isempty(CBIG_REPDATA_DIR))
    error('CBIG_REPDATA_DIR does not exist. Please check your config file.')
end
ref_results_dir = fullfile(CBIG_REPDATA_DIR, 'stable_projects', 'brain_parcellation', 'Yan2023_homotopic');

%% step 1: check if premultiplied matrix are identical
load(fullfile(ref_results_dir, 'GSP_fullset_full_corr_single.mat'), 'final_PMM', 'dim');
ref_dim = dim;
ref_final_PMM = final_PMM;

user_PMM = fullfile(output_dir, 'premultiplied_matrix_single.mat');
assert(logical(exist(user_PMM, 'file')), 'User failed to generate premultiplied matrix.');
load(user_PMM, 'final_PMM', 'dim');
user_dim = dim;
user_final_PMM = final_PMM;

pmm_diff = mean(abs(ref_final_PMM(:) - user_final_PMM(:)));

if(pmm_diff <= 1e-5 && ref_dim == user_dim)
    disp('Successfully generated the premultiplied matrix and the result matched the reference.');
else
    disp('FAIL!')
    disp('The mean difference between your results and reference results is ...');
    disp(['Premultiplied matrix difference: ' num2str(pmm_diff)]);
end

%% step 2: check if resultant parcellation is matched
load(fullfile(ref_results_dir,...
    '400parcels_C1.4e+05_K15_Wxyz1.5e+05_D50000_A1_iterations_100_seed_835.mat'), 'results');
ref_label = results.full_label;

user_label_dir = fullfile(output_dir, 'results',...
'400parcels_C1.4e+05_K15_Wxyz1.5e+05_D50000_A1_iterations_100_seed_835.mat');
assert(logical(exist(user_label_dir, 'file')), 'User failed to generate parcellation labels.');
load(user_label_dir, 'results');
user_label = results.full_label;

total_label_diff_count = sum(ref_label ~= user_label);
if(total_label_diff_count < 10) % there might be minor differences due to running environment
    disp('Successfully generated the level 400 parcellation and the result matched the reference.');
else
    disp('FAIL!')
    disp(['The total number of different labels between the user and reference labels is '...
        num2str(total_label_diff_count)]);
end
end
    