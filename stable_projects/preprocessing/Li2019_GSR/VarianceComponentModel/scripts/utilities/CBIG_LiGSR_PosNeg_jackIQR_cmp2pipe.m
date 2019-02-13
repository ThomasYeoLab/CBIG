function [IQR_pos, IQR_neg, med_pos, med_neg, mean_m2] = ...
    CBIG_LiGSR_PosNeg_jackIQR_cmp2pipe( m_del_d_mat1, m_del_d_mat2 )

% [IQR_pos, IQR_neg, med_pos, med_neg, mean_m2] = ...
%     CBIG_LiGSR_PosNeg_jackIQR_cmp2pipe( m_del_d_mat1, m_del_d_mat2 )
% 
% This function counts the number of behavioral measures where the entire
% interquartile range of the improved explained variance (preprocessing
% pipeline 1 versus pipeline 2) was above (or below) 0, and counts the
% number of behaviors where the median of the improved explained variance
% was above (or below) 0.
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
% Outputs:
%   - IQR_pos
%     The total number of traits where the whole interquartile range is above
%     zero, i.e. in at least 75% jackknife samples, pipeline 1 > pipeline 2.
% 
%   - IQR_neg
%     The total number of traits where the whole interquartile range is below
%     zero, i.e. in at least 75% jackknife samples, pipeline 2 < pipeline 1.
% 
%   - med_pos
%     The number of traits where the median value of the difference between
%     m_del_d_mat1 and m_del_d_mat2 across all jackknife samples is above 0
%     i.e. for each trait, check if median(m_del_d_mat1 and m_del_d_mat2) > 0.
% 
%   - med_neg
%     The number of traits where the median value of the difference between
%     m_del_d_mat1 and m_del_d_mat2 across all jackknife samples is below 0
%     i.e. for each trait, check if median(m_del_d_mat1 and m_del_d_mat2) < 0.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

dif = m_del_d_mat1 - m_del_d_mat2;
quart_25 = quantile(dif, 0.25, 1);
quart_75 = quantile(dif, 0.75, 1);
med = median(dif, 1);
mean_m2 = [mean(m_del_d_mat1(:)) mean(m_del_d_mat2(:))];

threshold = 1e-7;
IQR_pos = length(find(quart_25>-threshold));
IQR_neg = length(find(quart_75<threshold));
med_pos = length(find(med>0));
med_neg = length(find(med<0));

end

