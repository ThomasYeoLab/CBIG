function y = CBIG_MFMem_rfMRI_simBOLD_downsampling(x,bin)

%--------------------------------------------------------------------------
% y = CBIG_MFMem_rfMRI_simBOLD_downsampling(x,bin)
%
% Function for reducing data samples of x {n x samples} columns , i.e. bin=8, 
% put 8 points together, only count the first number of the 8 points. if the 
% residual is less than 8 points, count the first number of the residual.
%
% Input:
%     - x:      input data matrix {n x samples}
%     - bin:    put "bin" points in samples together, 
%
%  Output
%      - y:     reduced data matrix   
%
% 
%  Example:
%    y = CBIG_MFMem_rfMRI_simBOLD_downsampling(x,bin)
%    suppose:
%          - x = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18;
%                 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]                  ]
%          - bin = 4
%    y = CBIG_MFMem_rfMRI_simBOLD_downsampling(x,bin)
%     
%    1 2 3 4 -> 1;  5 6 7 8 -> 5; 9 10 11 12 -> 9; 13 14 15 16 -> 13;
%    17 18 -> 17,
%
%    y = [1 5 9 13 17; 1 5 9 13 17]
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%--------------------------------------------------------------------------
n = size(x,1);

if mod(size(x,2),bin) == 0
   T = size(x,2)/bin;
else
   T = fix(size(x,2)/bin)+1;
end

y = zeros(n,T);

for i = 1:T
    y(:,i) = x(:,bin*(i-1)+1);
end

end
