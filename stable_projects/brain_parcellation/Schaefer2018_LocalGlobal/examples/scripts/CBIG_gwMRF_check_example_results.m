function [if_both_hemi_match] = CBIG_gwMRF_check_example_results(output_folder, seed)
% This function checks if the generated example results are identical with
% the reference files provided for a given seed.
%
% Input:
%   - output_folder: The path under which the generated example results is saved.
%   - seed: indicate the index of the seed of whose example results will be compared.
%
% Output:
%   - if_both_hemi_match: logical variable indicates if the label files for both hemispheres match well.
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% check if replacing example results
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));

if(~ischar(seed))
    seed = num2str(seed);
end

example_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','brain_parcellation',...
'Schaefer2018_LocalGlobal','examples');

% make sure that the user-generated file exists
my_mat = dir(fullfile(output_folder,'clustering',['*_seed_', seed '.mat']));
if(isempty(my_mat))
    error('The input file %s is not found.', fullfile(output_folder,'clustering',['*_seed_', seed '.mat']));
else
    my_seed = load(fullfile(output_folder, 'clustering', my_mat.name));
end

ref_mat = dir(fullfile(example_dir,'example_results', 'clustering', ['*_seed_' seed '.mat']));
if(isempty(ref_mat))
    error('The example file %s is not found.',fullfile(example_dir,'example_results',...
        'clustering', ['*_seed_' seed '.mat']));
else
    ref_seed = load(fullfile(example_dir,'example_results','clustering',ref_mat.name));
end

if_both_hemi_match = true;

% check if the number of non medial wall voxels in lh ref labels equals the
% number of total overlapping voxels between the two lh label files.
ref_lh_labels = ref_seed.results.lh_label;
[~, ~, overlap_lh, ~] = CBIG_HungarianClusterMatch(ref_lh_labels, my_seed.results.lh_label, 1);
lh_ref_voxel_count = length(find(ref_lh_labels(:)~=0));

% check lh magnitude of overlapping
lh_diff = lh_ref_voxel_count - abs(overlap_lh);
if(lh_diff == 0)
    fprintf('Left hemisphere labels are the same.\n');
elseif(lh_diff/length(lh_ref_voxel_count) < 0.0005) % Less than 0.05% of total voxels are different.
    fprintf([num2str(lh_diff) ' left hemisphere are different.\n']);
    warning('The small difference might be caused by different running environments.');
else % Too many voxels are different
    fprintf([num2str(lh_diff) ' left hemisphere labels are different.\n']);
    if_both_hemi_match = false;
end

% now check for rh
ref_rh_labels = ref_seed.results.rh_label;
[~, ~, overlap_rh, ~] = CBIG_HungarianClusterMatch(ref_rh_labels, my_seed.results.rh_label, 1);
rh_ref_voxel_count = length(find(ref_rh_labels(:)~=0));

% check rh magnitude of overlapping
rh_diff = rh_ref_voxel_count - abs(overlap_rh);
if(rh_diff == 0)
    fprintf('Right hemisphere labels are the same.\n');
elseif(rh_diff/length(rh_ref_voxel_count) < 0.0005) % Less than 0.05% of total voxels are different.
    fprintf([num2str(rh_diff) ' right hemisphere are different.\n']);
    warning('The small difference might be caused by different running environments.');
else % Too many voxels are different
    fprintf([num2str(rh_diff) ' right hemisphere labels are different.\n']);
    if_both_hemi_match = false;
end

% check if we need to replace example results
if(replace_unittest_flag)
    disp('Replacing example reference results for Schaefer2018 LocalGlobal Parcellation...');
    copyfile(fullfile(output_folder, 'clustering', my_mat.name),...
     fullfile(example_dir,'example_results','clustering',ref_mat.name));
else
    % produce final parcellation visualization if both hemisphere match
    if(if_both_hemi_match)
        fprintf('Hurray! Example results for both hemispheres match.\n');
        CBIG_DrawSurfaceMaps(my_seed.results.lh_label, my_seed.results.rh_label, 'fsaverage6', 'inflated');
    end
end
end 
