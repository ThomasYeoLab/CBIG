function [perc_improv, m_jack, v_jack] = CBIG_LiGSR_del_d_jack_cmp2pipelines( ...
    m_del_d_mat1, m_del_d_mat2, num_families, d )

% [perc_improv, m_jack, v_jack] = CBIG_LiGSR_del_d_jack_cmp2pipelines( ...
%     m_del_d_mat1, m_del_d_mat2, num_families, num_samples, d )
% 
% This function computes three statistics of the mean explained variance
% difference between two preprocessing pipelines, averaged across all
% traits:
% (1) the percentage improvement of pipeline 1 relative to pipeline 2;
% (2) the mean of all jackknife samples;
% (3) the jackknife variance.
% 
% Inputs:
%   - m_del_d_mat1
%     A #JackknifeSamples x #traits matrix. Each element of the matrix
%     contains the explained variance estimated on each jackknife sample
%     for the first preprocessing pipeline.
% 
%   - m_del_d_mat2
%     A #JackknifeSamples x #traits matrix. Each element of the matrix
%     contains the explained variance estimated on each jackknife sample
%     for the second preprocessing pipeline.
% 
%   - num_families
%     A scalar, the number of families in the full subject list. If all
%     subjects are unrelated, then the number of families equals to the
%     number of subjects.
%  
%   - d
%     A scalar, the number of subjects that were removed for each jackknife
%     sample.
% 
% Outputs
%   - perc_improv
%     The percentage improvement of the mean explained behavioral variance
%     of pipeline 1 relative to pipeline 2.
% 
%   - m_jack
%     The jackknife mean of the explained variance difference between the
%     two pipelines (averaged across all traits).
% 
%   - v_jack
%     The jackknife variance of the explained variance difference between
%     the two pipelines (averaged across all traits).
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

m_delete_d = mean(m_del_d_mat1 - m_del_d_mat2, 2);     % num_samples x 1
num_samples = length(m_delete_d);
perc_improv = mean(m_delete_d) / mean(mean(m_del_d_mat2, 2), 1);

% jackknife mean and variance
m_jack = mean(m_delete_d, 1);            % 1 x 1
v_jack = (num_families - d) / (num_samples * d) * sum(bsxfun(@minus, m_delete_d, m_jack).^2, 1);   % 1 x 1


end

