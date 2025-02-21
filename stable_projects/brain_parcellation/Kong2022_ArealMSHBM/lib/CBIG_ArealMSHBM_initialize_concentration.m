function lambda = CBIG_ArealMSHBM_initialize_concentration(D)

% lambda = CBIG_ArealMSHBM_initialize_concentration(D)
%
% This function finds initialize point for concentration parameter with given feature dimension
%
% Input:
%   - D: feature dimension
%
% Output:
%   - lambda: initialize point for a given feature dimension D
%
% Example:
%   lambda = CBIG_ArealMSHBM_initialize_concentration(1482)
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


ini_flag = 0;
start_value = 50;
ini_loop = 0;
while (ini_flag == 0)
    ini_loop = ini_loop + 1;
    [inilambda,~,exitflag] = fzero(@(inputx) sign(inputx)*abs(besseli(D/2 - 1,inputx))-1e+10, ...
        start_value,optimset('Display','off'));
    if (exitflag == 1)
        lambda = inilambda;
        ini_flag = 1;
    else
        start_value = start_value + 50;
    end
    if(ini_loop == 1000)
        error('Can not initialize lambda automatically, please manually set it');
    end
end
fprintf('initialize lambda with %.2f\n',lambda);
        
end