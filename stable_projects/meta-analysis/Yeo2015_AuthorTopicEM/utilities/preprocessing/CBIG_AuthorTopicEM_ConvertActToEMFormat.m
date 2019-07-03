function formatted_act = CBIG_AuthorTopicEM_ConvertActToEMFormat(act)
% formatted_act = CBIG_AuthorTopicEM_ConvertActToEMFormat(act)
%
% Convert activation data to input format required for the
% Expectation-Maximization (EM) algorithm
%
% Input:
%   - act = activation foci across experiments
% Output:
%   - formatted_act  = E x 3 cell array of experiments' activation foci.
%     formatted_act{e, 1} = Ne x V sparse matrix where Ne is the number of unique activation foci
%                   in the given experiment, and V is the number of brain voxels.
%                   w{1} is the same as w{2} but has counts.
%     formatted_act{e, 2} = Ne x V sparse matrix where formatted_act{e, 2}(n, v) = 1 if the n-th activation focus in
%                   the experiment is the v-th voxel.
%     formatted_act{e, 3} = 1 x Nd vector, where formatted_act{n} is the number of times the n-th unique activation
%                 is reported in the experiment e.
%
% Example:
%   formatted_act = CBIG_AuthorTopicEM_ConvertActToEMFormat(act)
%   Convert activation data to input format for EM algorithm
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

formatted_act = cell(size(act, 1), 3);
V = size(act, 2);
for d = 1:length(formatted_act)
    unique_words = find(act(d, :) ~= 0);
    num_repeats = act(d, act(d, :) ~= 0);
    Nd = length(unique_words);
    
    formatted_act{d, 2} = sparse(transpose(1:Nd), unique_words(:), ones(Nd, 1), Nd, V);
    formatted_act{d, 1} = sparse(transpose(1:Nd), unique_words(:), num_repeats', Nd, V);
    formatted_act{d, 3} = num_repeats; 
end
