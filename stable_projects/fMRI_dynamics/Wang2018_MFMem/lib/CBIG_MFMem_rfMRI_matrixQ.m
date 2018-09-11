function q = CBIG_MFMem_rfMRI_matrixQ(i,n,T)

%------------------------------------------------------------------------
% q = CBIG_MFMem_rfMRI_matrixQ(i,n,T)
% 
% function for generation specitial diagnal matrix used in estimation
%
% Input:
%     - i   is the index to selected i-th trail
%     - n   is total availble data trails
%     - T   is number of samples of each data trail
%
% Output:
%     - q   is a matrix, q:{nT  x nT}
%
% Example:
% q = CBIG_EM_Q(2,2,3)
%                       trails  samples  selected
% q =[0 0 0 0 0 0               T1
%     0 0 0 0 0 0        n1     T2         No
%     0 0 0 0 0 0               T3
%     0 0 0 1 0 0               T1       
%     0 0 0 0 1 0        n2     T2         Yes
%     0 0 0 0 0 1]              T3        
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%-----------------------------------------------------------------------

    q_dig = zeros(1,n*T);
    q_dig(T*(i-1)+1:T*(i-1)+T) = 1;
    q = diag(q_dig);
end
