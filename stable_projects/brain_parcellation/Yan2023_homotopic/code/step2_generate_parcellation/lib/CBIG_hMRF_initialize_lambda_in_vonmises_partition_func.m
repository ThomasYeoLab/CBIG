function lambda = CBIG_hMRF_initialize_lambda_in_vonmises_partition_func(lambda_init, dimension_for_bessel_func)
% lambda = CBIG_hMRF_initialize_lambda_in_vonmises_partition_func(lambda_init, dimension_for_bessel_func)
%
% This function determines the initial value for lambda in the partition function of von Mises distribution.
% The initialization is determined by the algorithm if user passes in lambda_init == 0.
% Automatic initialization is recommended to avoid numerical overflow encountered at the bessel function.
% In this project, lambda corresponds to the concentration parameter in the vonmises function:
% 1) kappa of the normalization term in global likelihood term
% 2) tau of the normalization term in xyz likelihood term

% Input
%   - lambda_init: (double)
%     The initial lambda. Pass 0 if automatic initialization is desired.
%   - dimension_for_bessel_func: (double)
%     Dimensionality for the bessel function.
%
% Output
%   - lambda (double):
%     The initialized lambda.
%
% Example
%   - lambda = CBIG_hMRF_initialize_lambda_in_vonmises_partition_func(0, 3)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(lambda_init == 0)
    ini_flag = 0;
    if(dimension_for_bessel_func < 10)
        start_value = 50;
    elseif(dimension_for_bessel_func < 1000)
        start_value = 100;
    elseif(dimension_for_bessel_func < 4e5)
        start_value = 1e5;
    else
        start_value = 1e6;
    end
    ini_loop = 0;
    while (ini_flag == 0)
        ini_loop = ini_loop + 1;
        try
            [inilambda, ~, exitflag] = fzero(@(inputx)...
                sign(inputx)*abs(besseli(dimension_for_bessel_func,inputx))-1e+10,...
                 start_value,optimset('Display','off'));
        catch ME
            warning(['CBIG_hMRF_initialize_lambda_in_vonmises_partition_func cannot automatically initialize lambda'...
            ' for your input dimension. You can try another start_value by modifying line 25-33 '...
            'in CBIG_hMRF_initialize_lambda_in_vonmises_partition_func.']);
            rethrow(ME)
        end
        if (exitflag == 1)
            lambda = inilambda;
            ini_flag = 1;
        else
            start_value = start_value + 100;
        end
        if(ini_loop == 10000)
            error('Can not initialize lambda automatically, please manually set it');
        end
    end
else
    lambda = lambda_init;
end
end