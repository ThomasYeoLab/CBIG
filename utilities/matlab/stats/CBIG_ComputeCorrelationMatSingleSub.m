function CBIG_ComputeCorrelationMatSingleSub(src_mask, src_fMRI, targ_mask, targ_fMRI, output_file)

% Compute correlation matrix between source and target fMRI for a single subject.
%   CBIG_ComputeCorrelationMatSingleSub(src_mask, src_fMRI, targ_mask, targ_fMRI, output_file)
%   Input:
%       src_mask   : mask of source
%       src_fMRI   : fMRI data of source
%       targ_mask  : mask of target
%       targ_fMRI  : fMRI data of target
%       output_file: output file
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


s_mask = MRIread(src_mask);
s_fMRI = MRIread(src_fMRI);
s_series = zeros(sum(s_mask.vol(:) == 1), size(s_fMRI.vol, 4));

for t = 1:size(s_series, 2)
   tmp = squeeze(s_fMRI.vol(:, :, :, t)); 
   s_series(:, t) = tmp(s_mask.vol(:) == 1);
end
clear s_mask;
clear s_fMRI;
s_series = transpose(s_series);

t_mask = MRIread(targ_mask);
t_fMRI = MRIread(targ_fMRI);
t_series = zeros(sum(t_mask.vol(:) == 1), size(t_fMRI.vol, 4));
for t = 1:size(t_series, 2)
   tmp = squeeze(t_fMRI.vol(:, :, :, t)); 
   t_series(:, t) = tmp(t_mask.vol(:) == 1);
end
clear t_mask
clear t_fMRI
t_series = transpose(t_series);

% normalize series (note that series are now of dimensions: T x N)
s_series = s_series - repmat(mean(s_series, 1), size(s_series, 1), 1);
s_series = s_series./repmat(sqrt(sum(s_series.^2, 1)), size(s_series, 1), 1);

t_series = t_series - repmat(mean(t_series, 1), size(t_series, 1), 1);
t_series = t_series./repmat(sqrt(sum(t_series.^2, 1)), size(t_series, 1), 1);

corr_mat = zeros(size(s_series, 2), size(t_series, 2));
tic
for i = 1:size(s_series, 2)
  if(mod(i, 1000) == 0)
    disp(num2str(i)); toc
  end
    corr_mat(i, :) = sum(repmat(s_series(:, i), 1, size(t_series, 2)) .* t_series, 1);
end
toc

  save(output_file, 'corr_mat', '-v7.3');
exit
