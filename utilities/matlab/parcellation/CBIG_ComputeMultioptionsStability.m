function CBIG_ComputeMultioptionsStability(multi_mat_file, subset, exit_flag)

% CBIG_ComputeMultioptionsStability(multi_mat_file, subset, exit_flag)
%
% Wrapper function for computing the stability of clustering results in
% multiple option choices.
%
% "multi_mat_file" should have "lh_labels", 'rh_labels', and "likelihood"
% variables inside.
%
% The output stability file is located in the same folder as
% "multi_mat_file" with ".stability" suffix.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md




if(ischar(exit_flag))
   exit_flag = str2num(exit_flag); 
end

if(ischar(subset))
   subset = str2num(subset); 
end

load(multi_mat_file);
n = size(lh_labels, 1);
lh_index = find(lh_labels(:, 1) ~=0);
rh_index = find(rh_labels(:, 1) ~=0);

lh_labels = lh_labels(lh_index, :);
rh_labels = rh_labels(rh_index, :);
labels = [lh_labels; rh_labels];

% select subset
[likelihood_sorted, I] = sort(likelihood, 'descend');
labels = labels(:, I(1:subset));

stability = CBIG_ComputeClusterStability(labels);
lh_s = zeros(n, 1); lh_s(lh_index) = stability(1:length(lh_index));
rh_s = zeros(n, 1); rh_s(rh_index) = stability(length(lh_index)+1:end);

a1 = strfind(multi_mat_file, '.mat');
output_prefix = multi_mat_file(1:a1-1);
output_file = [output_prefix '.stability.' num2str(subset) '.mat'];
save(output_file, 'lh_s', 'rh_s');

if(exit_flag)
  exit
end