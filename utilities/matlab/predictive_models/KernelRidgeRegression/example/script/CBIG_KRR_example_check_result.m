function CBIG_KRR_example_check_result(user_output_path)

% CBIG_KRR_example_check_result(user_output_path)
% This fucntion compares the user-generated results with the reference
% results for each split (3 splits in total). 
% If the optimal accuracies differ too much from the allowable threshold,
% an error message will be generated.
% Input:
%  - user_output_path:
%    the absolute path of the output directory containing user-generated
%    results
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

for split = 1:3
    split = num2str(split);
    
    % load reference result
    ref_output_path = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'predictive_models', ...
        'KernelRidgeRegression', 'example', 'reference_output', split, ['final_result_CO_' split '.mat']);
    ref = load(ref_output_path);

    % load user-generated result
    user_output_final_result = fullfile(user_output_path, split, ['final_result_CO_' split '.mat']);
    load(user_output_final_result);

    % compare optimal accuracies
    diff = sum(abs(ref.optimal_acc - optimal_acc));
    if diff >= 1e-15
        fprintf('The optimal accuracies for split %s differs too much from the reference output.\n', split);
    else
        fprintf('The optimal accuracies for split %s agree with the reference output.\n', split)
    end
end

end

