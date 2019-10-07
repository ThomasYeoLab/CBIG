function avgCorr = CBIG_ASDf_indivCorr2avgCorr(lh2lh_corr, lh2rh_corr, ...
rh2rh_corr, lh2subcor_corr, rh2subcor_corr, subcor2subcor_corr)
% avgCorr = CBIG_ASDf_indivCorr2avgCorr(lh2lh_corr, lh2rh_corr, 
% rh2rh_corr, lh2subcor_corr, rh2subcor_corr, subcor2subcor_corr)
% 
% This function averages hemispheric/subcortical FC correlation matrices acorss
% subjects, and concatenates the hemispheric/subcortical correlation matrices together.
% 
% Input:
%     - lh2lh_corr:
%           MxMxN matrix, where M is the number of ROIs in the left/right
%           hemisphere, N is the number of subjects
%     - lh2rh_corr:
%           MxMxN matrix, where M is the number of ROIs in the left/right
%           hemisphere, N is the number of subjects
%     - rh2rh_corr:
%           MxMxN matrix, where M is the number of ROIs in the left/right
%           hemisphere, N is the number of subjects
%     - lh2subcor_corr:
%           MxPxN matrix, where M is the number of ROIs in the left/right
%           hemisphere, P is the number of subcortical ROIs, N is the number of subjects
%     - rh2subcor_corr:
%           MxPxN matrix, where M is the number of ROIs in the left/right
%           hemisphere, P is the number of subcortical ROIs, N is the number of subjects
%     - subcor2subcor_corr:
%           PxPxN matrix, where P is the number of subcortical ROIs, N is the number of subjects
%
% Output:
%     - avgCorr:
%           (2*M+P)x(2*M+P) matrix. Left hemisphere first, followed by right hemisphere. 
%           Brain networks order is not re-arranged.
% 
% Example:
%       avgCorr = CBIG_ASDf_indivCorr2avgCorr(lh2lh_corrmat, lh2rh_corrmat, 
%       rh2rh_corrmat, lh2subcor_corrmat, rh2subcor_corrmat, subcor2subcor_corrmat)
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

for idx = 1:size(lh2lh_corr,3)
    lh2lh_corr(:,:,idx) = CBIG_StableAtanh(lh2lh_corr(:,:,idx));
    lh2rh_corr(:,:,idx) = CBIG_StableAtanh(lh2rh_corr(:,:,idx));
    rh2rh_corr(:,:,idx) = CBIG_StableAtanh(rh2rh_corr(:,:,idx));
    lh2subcor_corr(:,:,idx) = CBIG_StableAtanh(lh2subcor_corr(:,:,idx));
    rh2subcor_corr(:,:,idx) = CBIG_StableAtanh(rh2subcor_corr(:,:,idx));
    subcor2subcor_corr(:,:,idx) = CBIG_StableAtanh(subcor2subcor_corr(:,:,idx));
end

lh2lh_avgCorr = tanh(mean(lh2lh_corr,3));
lh2rh_avgCorr = tanh(mean(lh2rh_corr,3));
rh2rh_avgCorr = tanh(mean(rh2rh_corr,3));
lh2subcor_avgCorr = tanh(mean(lh2subcor_corr,3));
rh2subcor_avgCorr = tanh(mean(rh2subcor_corr,3));
subcor2subcor_avgCorr = tanh(mean(subcor2subcor_corr,3));

avgCorr = [lh2lh_avgCorr lh2rh_avgCorr lh2subcor_avgCorr; lh2rh_avgCorr' rh2rh_avgCorr rh2subcor_avgCorr; ...
lh2subcor_avgCorr' rh2subcor_avgCorr' subcor2subcor_avgCorr];
avgCorr = (avgCorr+avgCorr')/2;
