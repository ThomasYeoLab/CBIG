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
ref = matfile(fullfile(ref_results_dir, 'GSP_fullset_full_corr_single.mat'));
ref_dim = ref.dim;

user_PMM_dir = fullfile(output_dir, 'premultiplied_matrix_single.mat');
assert(logical(exist(user_PMM_dir, 'file')), 'User failed to generate premultiplied matrix.');

user = matfile(user_PMM_dir);
user_dim = user.dim;

pmm_diff = 0;
% performing partial check (1 row from every 100 rows), since it would take
% hours to check the entire PMM.
for i = 1:50:length(ref.final_PMM(:,1))
    ref_row = ref.final_PMM(i,:);
    user_row = user.final_PMM(i,:);
    pmm_diff = pmm_diff + mean(abs(ref_row - user_row));
end

if(pmm_diff <= 1e-6 && ref_dim == user_dim)
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
    