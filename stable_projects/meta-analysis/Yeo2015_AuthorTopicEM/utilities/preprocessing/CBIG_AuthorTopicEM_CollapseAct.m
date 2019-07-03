function [act, paradigm_by_exp, num_studies_per_paradigm] = CBIG_AuthorTopicEM_CollapseAct(act, paradigm_by_exp)
% [act, paradigm_by_exp, num_studies_per_paradigm] = CBIG_AuthorTopicEM_CollapseAct(act, paradigm_by_exp)
%
% Summarize input data by unique task paradigms to reduce storage size for the variables
%
% Input:
%   - act = E x V array of activation, where E is the number of experiments and V is the
%           number of brain voxels.
%   - paradigm_by_exp = E x T binary array, where E is the number of
%           experiments, and T is the number of task paradigms.
% Output:
%   - act = new array of activation with number of rows equals the number of unique
%           task paradigms
%   - paradigm_by_exp = E x T* binary array, where T* is the number of unique
%           paradigms
%   - num_studies_per_paradigm = T* x 1 array, with the t-th element contains the
%     number of experiments that recruits the unique task paradigm t
%
% Example:
%   [act, paradigm_by_exp, num_studies_per_paradigm] = CBIG_AuthorTopicEM_CollapseAct(...
%     act, paradigm_by_exp)
%   Collapse the original data in act and paradigm_by_exp
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  num_paradigms = size(paradigm_by_exp, 1);
  binary_vec = [2.^((num_paradigms-1):-1:0)]';
  paradigms  = sum(bsxfun(@times, paradigm_by_exp, binary_vec), 1);
  
  [unique_paradigm_comb, I] = unique(paradigms, 'first');
  paradigm_by_exp = paradigm_by_exp(:, I);
  
  new_act = zeros(length(unique_paradigm_comb), size(act, 2));
  for i = 1:length(unique_paradigm_comb)
      num_studies_per_paradigm(i) = sum(paradigms == unique_paradigm_comb(i));
      new_act(i, :) = sum(act(paradigms == unique_paradigm_comb(i), :), 1);
  end
  act = new_act;

