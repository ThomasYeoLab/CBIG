function null_hist = CBIG_GenerateALENullDist_EickhoffUnionAnalytical(ale_mat, bin_width)

% null_hist = CBIG_GenerateALENullDist_EickhoffUnionAnalytical(ale_mat, bin_width)
%
% Generate a histogram of ALE values from individual experiments under a null hypothesis of spatial independence
% FORMAT null_hist = CBIG_GenerateALENullDist_EickhoffUnionAnalytical(ale_mat, [bin_width])
%
% ale_mat = # studies x # relevant voxels matrix of activation probabilities
% bin_width = bin width of the histograms, set to 0.00001 if no argument was supplied 
% 
% null_hist = a histogram of ALE values, in which:
%         0 <= all values lie within null_hist(1) < bin_width
% bin_width <= all values lie within null_hist(2) < 2*bin_width
%               ...
% 
% ____________________________________________________________________________________________________
%
% CBIG_GenerateALENullDist_EickhoffUnionAnalytical returns a histogram of ALE values under a null-hypothesis of
% random, spatially independent brain locations as described in Eickhoff et. al, 2012.
% ____________________________________________________________________________________________________
%
% Refs:
%
% Eickhoff et al., NeuroImage, 2012. Activation likelihood estimation meta-analysis revisited
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if(nargin < 2)
      bin_width = 0.00001;
  end
  num_bins  = round(1/bin_width);

  % convert ale_mat to histogram
  hist_mat = zeros(size(ale_mat, 1), num_bins);
  for i = 1:size(ale_mat, 1)
     hist_mat(i, :) = ConvertProb2Bins(ale_mat(i, :), num_bins);
  end

  % convolve!
  null_hist = hist_mat(1, :);
  for i = 2:size(hist_mat, 1)
     null_hist = ConvolveALEHistogram(null_hist, hist_mat(i, :)); 
  end


% convert a vector of probabilities into a histogram of num_bins bins
function hist_vec = ConvertProb2Bins(prob_vec, num_bins)
  prob_vec = floor(prob_vec*num_bins)+1;
  [xi, m]  = knt2brk(prob_vec);

  hist_vec     = zeros(1, num_bins);
  hist_vec(xi) = m;
  hist_vec = hist_vec/sum(hist_vec);


% convolve two ALE histograms
function new_hist = ConvolveALEHistogram(hist_vec1, hist_vec2)
  index1 = find(hist_vec1 > 0);
  index2 = find(hist_vec2 > 0);

  num_bins = length(hist_vec2);
  new_hist = zeros(1, num_bins);
  for i = index1
      
      old_val1 = (i - 1)*1/num_bins;
      for j = index2
          
          old_val2 = (j - 1)*1/num_bins;
          new_val  = 1 - (1 - old_val1)*(1 - old_val2);
          
          new_hist(floor(new_val*num_bins)+1) = new_hist(floor(new_val*num_bins)+1) + hist_vec1(i)*hist_vec2(j);
      end
  end

  new_hist = new_hist/sum(new_hist);
